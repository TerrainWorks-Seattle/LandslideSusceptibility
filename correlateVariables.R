source("./shared/createInitiationRange.R")
source("./shared/createNoninitiationRaster.R")
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


# The region where non-initiation buffers can be sampled from: regions that 
# meet initiation conditions but recorded no landslides
noninitRaster <- createNoninitiationRaster(
  varsRaster,
  initRange,
  initBuffers
)

# Generate non-initiation buffers
noninitBuffersCount <- ceiling(length(initPoints) * noninitRatio)

noninitBuffers <- generateNoninitiationBuffers(
  noninitBuffersCount,
  noninitRaster,
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
