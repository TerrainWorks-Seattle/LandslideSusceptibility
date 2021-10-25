refRasterFile          <- "~/Work/Data/Scottsburg/elev_scottsburg.flt"
varRasterFiles         <- list(
  "~/Work/Data/Scottsburg/elev_scottsburg.flt"
)   
initPointsFile         <- "~/Work/Data/Scottsburg/Scottsburg_Upslope.shp"
noninitRatio           <- 1.5
bufferRadius           <- 20
bufferExtractionMethod <- "center cell"
initRangeExpansion     <- 0
iterations             <- 20
testingProportion      <- 10
outputDir              <- "~/Work/Data/Scottsburg"

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

terra::plot(analysisRegionMask)
terra::polys(initBuffers, col = "blue")
terra::polys(noninitBuffers, col = "red")

# Create dataset -------------------------------------------------------------

# By default, extract all values from initiation and non-initiation buffers
initValues <- terra::extract(varsRaster, initBuffers, xy = TRUE)
noninitValues <- terra::extract(varsRaster, noninitBuffers, xy = TRUE)

extractionMethod <- "center cell"

# Subset if a different buffer extraction method was requested
if (extractionMethod == "center cell") {
  
  # Find center points of buffers
  initCenters <- terra::centroids(initBuffers)
  noninitCenters <- terra::centroids(noninitBuffers)
  
  # Get center point coordinates
  initCoords <- terra::geom(initCenters)[,c("x","y")]
  noninitCoords <- terra::geom(noninitCenters)[,c("x","y")]
  
  # Extract values from cells containing center points
  initValues <- terra::extract(varsRaster, initCoords, xy = TRUE)
  noninitValues <- terra::extract(varsRaster, noninitCoords, xy = TRUE)
  
} else if (extractionMethod == "max gradient cell") {
  
  initValues <- aggregateBufferValues(initValues, "grad", max)
  noninitValues <- aggregateBufferValues(noninitValues, "grad", max)
  
} else if (extractionMethod == "max plan cell") {
  
  initValues <- aggregateBufferValues(initValues, "plan", max)
  noninitValues <- aggregateBufferValues(noninitValues, "plan", max)
  
}

# Assign a classification value to each entry
initValues$class <- rep("initiation", nrow(initValues))
noninitValues$class <- rep("non-initiation", nrow(noninitValues))

# Combine initiation and non-initiation entries into a single dataset
dataset <- rbind(initValues, noninitValues)

# Remove the "ID" column
dataset$ID <- NULL

# TODO: Remove cells that fall within initiation AND non-initiation buffers

# Factor the classification variable values
dataset$class <- factor(dataset$class)

# Filter out entries with NA values
dataset <- na.omit(dataset)



task <- mlr3spatiotempcv::TaskClassifST$new(
  "landslides",
  backend = mlr3::as_data_backend(dataset),
  target = "class",
  positive = "initiation",
  extra_args = list(coordinate_names = c("x", "y"), crs = terra::crs(refRaster, proj = TRUE))
)

# Define resampling method for repeated k-fold spatial cross-validation
resampling <- mlr3::rsmp("repeated_spcv_coords", folds = 5, repeats = 10)

mlr3spatiotempcv::autoplot(resampling, task, fold_id = c(1:5), repeats_id = 1)

# Define a random forest learner
learner <- mlr3::lrn("classif.randomForest", predict_type = "prob")

# Perform resampling method with the random forest learner and the landslides task
result <- mlr3::resample(
  task = task,
  learner = learner,
  resampling = resampling
)

# Display ROC plot for model 1 predictions
mlr3viz::autoplot(result$predictions()[[1]], type = "roc")

# Rows used for testing model 1
result$resampling$test_set(1)

# Predictions for model 1
result$predictions()[[1]]

# Confusion matrix for model 1 predictions
result$predictions()[[1]]$confusion

# Calculate average classification error
result$aggregate(measures = mlr3::msr("classif.ce"))

# Calculate average area under ROC curve
result$aggregate(measures = mlr3::msr("classif.auc"))

# Calculate average accuracy (proportion of correct predictions out of all predictions)
result$aggregate(measures = mlr3::msr("classif.acc"))
