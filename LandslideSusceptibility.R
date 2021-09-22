referenceRasterFile <- "C:/Work/netmapdata/Scottsburg/elev_scottsburg.flt"
variableRasterFiles <- list(
  "C:/Work/netmapdata/Scottsburg/grad_30.tif",
  "C:/Work/netmapdata/Scottsburg/plan_15.tif"
)
initiationPointsFile <- "C:/Work/netmapdata/Scottsburg/Scottsburg_Upslope.shp"
bufferRadius <- 20
initiationLimitScaler <- 1.1

# Load rasters -----------------------------------------------------------------

# Load reference raster
referenceRaster <- terra::rast(referenceRasterFile)

# Load input rasters
inputRasters <- list(
  terra::rast(variableRasterFiles[[1]]),
  terra::rast(variableRasterFiles[[2]])
)

# Align input rasters
inputRasters <- TerrainWorksUtils::alignRasters(referenceRaster, inputRasters)

# Combine all input variable rasters into a single multi-layer raster
inputVarsRaster <- c(inputRasters[[1]], inputRasters[[2]])

# Create initiation regions ----------------------------------------------------

# Load initiation points
initiationPoints <- terra::vect(initiationPointsFile)

# Create buffer regions around initiation points
initiationPolys <- terra::buffer(initiationPoints, width = bufferRadius)

# Calculate variables' initiation range ----------------------------------------

# Extract input raster values within initiation regions
initiationPoints <- terra::project(initiationPoints, referenceRaster)
initiationValues <- terra::extract(inputVarsRaster, initiationPolys)
initiationValues <- initiationValues[-1]

# Define matrix to hold initiation limits
initiationRange <- matrix(rep(NA, 4), nrow = 2)
colnames(initiationRange) <- names(inputVarsRaster)
rownames(initiationRange) <- c("min", "max")

# Populate matrix with limit values
for (varName in names(inputVarsRaster)) {
  initiationRange["min", varName] <- min(initiationValues[[varName]], na.rm = TRUE)
  initiationRange["max", varName] <- max(initiationValues[[varName]], na.rm = TRUE)
}

# Scale limit values scale to broaden the limit
initiationRange <- initiationRange * initiationLimitScaler

# Create initiation mask -------------------------------------------------------

# Get all input variable values
inputVarValues <- terra::values(inputVarsRaster)

for (varName in names(inputVarsRaster)) {
  varRaster <- inputVarsRaster[[varName]]
  
  # Get variable value limits
  minInitiationValue <- initiationRange["min", varName]
  maxInitiationValue <- initiationRange["max", varName]
  
  # NA-out cells with values outside variable initiation range
  varRaster <- terra::app(varRaster, fun = function(x) {
    ifelse(x < minInitiationValue | x > maxInitiationValue, NA, x)
  })
  
  # Update the raster in the input raster stack
  inputVarsRaster[[varName]] <- varRaster
}

# NA-out all cells with values outside of initiation range
initiationMask <- terra::app(inputVarsRaster, fun = "prod")

# NA-out all cells inside initiation regions
initiationCells <- terra::extract(inputVarsRaster, initiationPolys, cells = TRUE)[["cell"]]
initiationMask[initiationCells] <- NA

# Generate non-initiation regions ----------------------------------------------

# Sample points in areas that fit initiation parameters but recorded no landslides
noninitiationPoints <- terra::spatSample(
  initiationMask, size = 10, na.rm = TRUE, as.points = TRUE
)

# Create buffer regions around non-initiation points
noninitiationPolys <- terra::buffer(noninitiationPoints, width = bufferRadius)

# Create training dataset ------------------------------------------------------

noninitiationValues <- terra::extract(inputVarsRaster, noninitiationPolys)
noninitiationValues <- noninitiationValues[-1]

initiationValues$class <- rep("initiation", nrow(initiationValues))
noninitiationValues$class <- rep("non-initiation", nrow(noninitiationValues))

# Combine initiation and non-initiation values into a single dataset
trainingData <- rbind(initiationValues, noninitiationValues)

trainingData <- na.omit(trainingData)

trainingData$class <- factor(trainingData$class)

# Create random forest model ---------------------------------------------------

rfModel <- randomForest::randomForest(
  formula = class ~ .,
  data = trainingData,
  ntree = 200
)

# Generate probability raster --------------------------------------------------

probRaster <- terra::predict(
  inputVarsRaster,
  rfModel,
  na.rm = TRUE,
  type = "prob"
)

