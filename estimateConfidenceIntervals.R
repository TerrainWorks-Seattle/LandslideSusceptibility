

refRasterFile          <- "E:/NetmapData/Scottsburg/elev_scottsburg.flt"
varRasterFiles         <- list(
                            "E:/NetmapData/Scottsburg/grad_30.tif",
                            "E:/NetmapData/Scottsburg/plan_30.tif",
                            "E:/NetmapData/Scottsburg/pca_scott.flt"
                          )   
initPointsFile         <- "E:/Netmapdata/Scottsburg/Scottsburg_Upslope.shp"
noninitRatio           <- 1.5
bufferRadius           <- 20
bufferExtractionMethod <- "center cell"
initRangeExpansion     <- 0
iterations             <- 20
testingProportion      <- 10
outputDir              <- "E:/NetmapData"

# Load helper functions ------------------------------------------------------

source("./shared/createInitiationRange.R")
source("./shared/createAnalysisRegionMask.R")
source("./shared/generateNoninitiationBuffers.R")
source("./shared/extractBufferValues.R")
source("./shared/createProportionRaster.R")

# Validate parameters --------------------------------------------------------

# Validate reference raster file
if (!file.exists(refRasterFile))
  stop(paste0("Reference raster file not found: '", refRasterFile, "'."))

# Validate variable raster files
lapply(varRasterFiles, function(file) {
  if (!file.exists(file))
    stop(paste0("Variable raster file not found: '", file, "'."))
})

# Validate initiation points file
if (!file.exists(initPointsFile))
  stop(paste0("Initiation points file not found: '", initPointsFile, "'."))

# Validate non-initiation ratio
if (noninitRatio <= 0)
  stop("Non-initiation points ratio cannot be <=0.")

# Validate buffer radius
if (bufferRadius < 0)
  stop("Buffer radius cannot be negative.")

# Validate buffer extraction method
if (!(bufferExtractionMethod %in% c("all cells", "center cell", "max gradient 
      cell", "max plan cell")))
  stop("Buffer extraction method must be one of 'all cells', 'center cell', 
       'max gradient cell', or 'max plan cell'.")

# Validate number of iterations
if (iterations < 1)
  stop("Iterations cannot be fewer than 1.")

# Validate testing proportion
if (testingProportion <= 0 || testingProportion >= 100)
  stop("Testing proportion must be between 0% and 100%.")

# Validate output directory
if (!file.exists(outputDir))
  stop("Output directory not found: '", outputDir, "'.")

# Load rasters ---------------------------------------------------------------

# Load reference raster
refRaster <- terra::rast(refRasterFile)

# Load explanatory variable rasters
varRasterList <- lapply(varRasterFiles, function(file) terra::rast(file))

# Align variable rasters
varRasterList <- TerrainWorksUtils::alignRasters(refRaster, varRasterList)

# Combine variable rasters into a single multi-layer raster
varsRaster <- terra::rast(varRasterList)

# Create initiation buffers --------------------------------------------------

# Load initiation points
initPoints <- terra::vect(initPointsFile)
initPoints <- terra::project(initPoints, refRaster)

# Create a buffer around each initiation point
initBuffers <- if (bufferRadius > 0) {
  terra::buffer(initPoints, width = bufferRadius)
} else {
  initPoints
}

# Generate non-initiation buffers --------------------------------------------

# Calculate initiation range for each variable
initRange <- createInitiationRange(
  varsRaster,
  initBuffers,
  initRangeExpansion
)

# Identify cells in the study region that have variable values within their 
# initiation ranges
analysisRegionMask <- createAnalysisRegionMask(
  varsRaster,
  initRange
)

# Define the region where non-initiation points can be generated

# NOTE: initiation buffers and non-initiation buffers must not overlap. This can
# be avoided by making sure each non-initiation point is generated at least 2 
# buffer-radius-lengths away from any initiation point 

# Double the size of the initiation buffers
expInitBuffers <- terra::buffer(initPoints, width = bufferRadius * 2)

# Remove expanded initiation buffers from the viable non-initiation region
noninitRegion <- terra::copy(analysisRegionMask)
initCellIndices <- terra::extract(noninitRegion, expInitBuffers, cells = TRUE)$cell
noninitRegion[initCellIndices] <- NA

# Determine how many non-initiation buffers to generate
noninitBuffersCount <- ceiling(length(initPoints) * noninitRatio)

# Generate non-initiation buffers
noninitBuffers <- generateNoninitiationBuffers(
  noninitBuffersCount,
  noninitRegion,
  bufferRadius
)

# Plot variable distributions ------------------------------------------------

# Get requested values from all initiation and non-initiation buffers
allBuffersData <- extractBufferValues(
  varsRaster,
  initBuffers,
  noninitBuffers,
  bufferExtractionMethod
)

# Create training & testing non-initiation buffer sets -----------------------

# Create a static set of training and testing non-initiation buffers to use every iteration
testingNoninitCount <- floor(length(noninitBuffers) * (testingProportion / 100))
testingNoninitIndices <- sample(seq_along(noninitBuffers), size = testingNoninitCount)
trainingNoninitIndices <- setdiff(seq_along(noninitBuffers), testingNoninitIndices)

trainingNoninitBuffers <- noninitBuffers[trainingNoninitIndices]
testingNoninitBuffers <- noninitBuffers[testingNoninitIndices]

# Perform cross-validation ---------------------------------------------------

probRasterList <- list()

# Calculate how many initiation buffers should be used for testing per iteration
testingInitBuffersCount <- floor(length(initBuffers) * (testingProportion / 100))

# Create a version of the variables raster that only keeps cell values within 
# the analysis region
analysisRegionVarsRaster <- terra::mask(varsRaster, analysisRegionMask)

for (i in seq_len(iterations)) {
  ## Create initiation buffer sets -------------------------------------------
  
  # Create testing initiation set
  testingInitIndices <- sample(seq_along(initBuffers), size = testingInitBuffersCount)
  testingInitBuffers <- initBuffers[testingInitIndices]
  
  # Create training initiation set
  trainingInitIndices <- setdiff(seq_along(initBuffers), testingInitIndices)
  trainingInitBuffers <- initBuffers[trainingInitIndices]
  
  ## Create model ------------------------------------------------------------
  
  # Create training dataset
  trainingData <- extractBufferValues(
    varsRaster,
    trainingInitBuffers,
    trainingNoninitBuffers,
    bufferExtractionMethod
  )
  
  # Train a new random forest model
  rfModel <- randomForest::randomForest(
    formula = class ~ .,
    data = trainingData
  )
  
  ## Generate probability raster ---------------------------------------------
  
  # Predict landslide probability values for all cells within the analysis
  # region
  probRaster <- terra::predict(
    analysisRegionVarsRaster,
    rfModel,
    na.rm = TRUE,
    type = "prob"
  )[["initiation"]]
  
  probRasterList[[i]] <- probRaster
  
}

if (length(probRasterList) > 1) {
  probsRaster <- terra::rast(probRasterList)
  
  firstPointCellIndex <- terra::extract(refRaster, initPoints[1], cell = TRUE)$cell
  cellProbs <- as.numeric(probsRaster[firstPointCellIndex])
  hist(cellProbs)
  
  probRasterAvg <- terra::app(probsRaster, fun = "mean")
  probRasterStd <- terra::app(probsRaster, fun = "std")
  
  probRasterLow <- probRasterAvg - probRasterStd
  probRasterHigh <- probRasterAvg + probRasterStd
  
  names(probRasterLow) <- "low"
  names(probRasterHigh) <- "high"
  
  propRasterLow <- createProportionRaster(probRasterLow)
  propRasterHigh <- createProportionRaster(probRasterHigh)
  
  terra::hist(probRasterStd, maxcell = 100000)
  terra::boxplot(probRasterStd, maxcell = 100000, horizontal = TRUE)  
}

filtCellIndices <- na.omit(which(terra::values(propRasterHigh)[,1] < 0.5))
propRasterHigh[filtCellIndices] <- NA
