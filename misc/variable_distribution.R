# Draws a boxplot for each input variable which compares the distribution of its 
# initiation and non-initiation values 

refRasterFile = "E:/NetmapData/Scottsburg/elev_scottsburg.flt"
varRasterFiles = list(
  "E:/NetmapData/Scottsburg/grad_30.tif",
  "E:/NetmapData/Scottsburg/plan_30.tif",
  "E:/NetmapData/Scottsburg/pca_scott.flt"
)
initPointsFile = "E:/NetmapData/Scottsburg/Scottsburg_Upslope.shp"
noninitRatio = 1
bufferRadius = 20
bufferExtractionMethod = "max gradient cell"
initRangeExpansion = 0

# Install dependencies ---------------------------------------------------------

dependencies <- c(
  "mlr3",
  "mlr3spatiotempcv",
  "terra"
)

install.packages(setdiff(dependencies, rownames(installed.packages())))

# Load helper functions --------------------------------------------------------

source("./helper/createInitiationRange.R")
source("./helper/createAnalysisRegionMask.R")
source("./helper/generateNoninitiationBuffers.R")
source("./helper/extractBufferValues.R")

# Load rasters -----------------------------------------------------------------

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

# Generate non-initiation buffers ----------------------------------------------

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

# Plot variable distributions --------------------------------------------------

# Get requested values from all initiation and non-initiation buffers
landslideData <- extractBufferValues(
  varsRaster,
  initBuffers,
  noninitBuffers,
  bufferExtractionMethod
)

# For each explanatory variable, draw a box plot of all the buffer values
for (varName in names(varsRaster)) {
  
  # get class values
  varInitValues <- landslideData[landslideData$class == "initiation", varName]
  varNoninitValues <- landslideData[landslideData$class == "non-initiation", varName]
  
  # Calculate class averages
  varInitAvg <- mean(varInitValues, na.rm = TRUE)
  varNoninitAvg <- mean(varNoninitValues, na.rm = TRUE)
  
  # Draw plot
  boxplot(
    varInitValues, varNoninitValues,
    main = paste0(varName, " distribution"),
    names = c("Initiation", "Non-initiation"),
    col = c("green", "red")
  )
  points(
    x = 1:2,
    y = c(varInitAvg, varNoninitAvg),
    pch = 16,
    cex = 2,
    col = "black"
  )
  
}