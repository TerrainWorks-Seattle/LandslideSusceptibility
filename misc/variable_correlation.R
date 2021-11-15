source("./shared/createInitiationRange.R")
source("./shared/createAnalysisRegionMask.R")
source("./shared/generateNoninitiationBuffers.R")
source("./shared/extractBufferValues.R")

refRasterFile = "E:/NetmapData/Scottsburg/elev_scottsburg.flt"
varRasterFiles = list(
  "E:/NetmapData/Scottsburg/grad_30.tif",
  "E:/NetmapData/Scottsburg/plan_30.tif",
  "E:/NetmapData/Scottsburg/pca_scott.flt"
)
initPointsFile = "E:/NetmapData/Scottsburg/Scottsburg_Upslope.shp"
noninitRatio = 2
bufferRadius = 20
bufferExtractionMethod = "center cell"
initRangeExpansion = 10

# -----------------------------------------------------------------------------

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

# Plot variable correlations ---------------------------------------------------

# Get data from all buffers
allBuffersData <- extractBufferValues(
  varsRaster,
  initBuffers,
  noninitBuffers,
  bufferExtractionMethod
)

coordsCols <- names(landslideData) %in% c("x", "y")  
allBuffersData <- allBuffersData[,!coordsCols]

# Plot each variable against the others
pairs(
  allBuffersData[1:3],
  col = ifelse(allBuffersData$class == "initiation", rgb(0,0,1,0.5), rgb(1,0,0,0.5)),
  pch = 16,
  cex = 0.8,
  lower.panel = NULL,
  main = paste0("Correlation for '", bufferExtractionMethod, "' of each buffer")
)
par(xpd = TRUE)
legend(
  "bottomleft",
  legend = c("Initiation", "Non-initiation"),
  col = c("blue", "red"),
  pch = 16
)
