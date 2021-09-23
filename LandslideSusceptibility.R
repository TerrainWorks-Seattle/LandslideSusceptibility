tool_exec <- function(in_params, out_params) {
  
  # Set parameters -------------------------------------------------------------
  
  # Input
  referenceRasterFile <- in_params[[1]]
  varRasterFiles <- in_params[[2]]
  initiationPointsFile <- in_params[[3]]
  bufferRadius <- in_params[[4]]
  initiationLimitScaler <- in_params[[5]]
  noninitiationPointsCount <- in_params[[6]]
  iterations <- in_params[[7]]
  outputDir <- in_params[[8]]
  
  # Set up logging -------------------------------------------------------------
  
  logFilename <- "landslide_susceptibility_log.txt"
  file.create(paste0(outputDir, "/", logFilename))
  
  logObj <- function(obj) {
    capture.output(obj, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  logMsg <- function(msg) {
    cat(msg, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  logMsg("Running with input parameters:\n")
  logMsg(paste0("  Reference raster: ", referenceRasterFile, "\n"))
  logMsg("  Explanatory variable rasters:\n")
  for (i in seq_along(varRasterFiles)) { logMsg(paste0("  [", i, "] ", varRasterFiles[[i]], "\n")) }
  logMsg(paste0("  Initiation points: ", initiationPointsFile, "\n"))
  logMsg(paste0("  Buffer radius: ", bufferRadius, "\n"))
  logMsg(paste0("  Initiation limit scaler: ", initiationLimitScaler, "\n"))
  logMsg(paste0("  Non-initiation points count: ", noninitiationPointsCount, "\n"))
  logMsg(paste0("  Iterations: ", iterations, "\n"))
  logMsg(paste0("  Output directory: ", outputDir, "\n"))
  logMsg("\n")
  
  # Load rasters ---------------------------------------------------------------
  
  # Load reference raster
  referenceRaster <- terra::rast(referenceRasterFile)
  logMsg(paste0("Loaded reference raster: ", referenceRasterFile, "\n"))
  
  # Load explanatory variable rasters
  varRasterList <- lapply(varRasterFiles, function(file) terra::rast(file))
  
  # Align variable rasters
  varRasterList <- TerrainWorksUtils::alignRasters(referenceRaster, varRasterList)
  
  # Combine variable rasters into a single multi-layer raster
  varsRaster <- terra::rast(varRasterList)
  
  # Create initiation regions --------------------------------------------------
  
  # Load initiation points
  initiationPoints <- terra::vect(initiationPointsFile)
  logMsg(paste0("Loaded initiation points: ", initiationPointsFile, "\n"))
  logMsg("\n")
  
  # Create buffer regions around initiation points
  initiationPolys <- terra::buffer(initiationPoints, width = bufferRadius)
  
  # Calculate variables' initiation range --------------------------------------
  
  # Extract variable values within initiation regions
  initiationPoints <- terra::project(initiationPoints, varsRaster)
  initiationValues <- terra::extract(varsRaster, initiationPolys)
  initiationValues <- initiationValues[-1] # Remove "ID" column
  
  # Define matrix to hold each variable's initiation value range
  initiationRange <- matrix(rep(NA, 4), nrow = 2)
  colnames(initiationRange) <- c("min", "max")
  rownames(initiationRange) <- names(varsRaster)
  
  # Populate matrix with range limits
  for (varName in names(varsRaster)) {
    initiationRange[varName, "min"] <- min(initiationValues[[varName]], na.rm = TRUE)
    initiationRange[varName, "max"] <- max(initiationValues[[varName]], na.rm = TRUE)
  }
  
  # Expand limit values scale to broaden the limit
  initiationRange <- initiationRange * initiationLimitScaler
  
  logMsg("Landslide initiation ranges:\n")
  logObj(initiationRange)
  
  # Create initiation mask -----------------------------------------------------
  
  # Define a raster which will NA-out any cells with variable values outside of 
  # their initiation ranges
  initiationRaster <- terra::copy(varsRaster)
  
  for (varName in names(initiationRaster)) {
    varRaster <- initiationRaster[[varName]]
    
    # Get variable value limits
    minInitiationValue <- initiationRange[varName, "min"]
    maxInitiationValue <- initiationRange[varName, "max"]
    
    # NA-out cells with values outside variable initiation range
    varInitiationRaster <- terra::app(varRaster, fun = function(x) {
      ifelse(x < minInitiationValue | x > maxInitiationValue, NA, x)
    })
    
    # Update the raster in the input raster stack
    initiationRaster[[varName]] <- varInitiationRaster
  }
  
  # NA-out all cells with variable values outside of their initiation range
  initiationRaster <- terra::app(initiationRaster, fun = "prod")
  
  # NA-out all cells inside initiation regions
  initiationCells <- terra::extract(initiationRaster, initiationPolys, cells = TRUE)[["cell"]]
  initiationRaster[initiationCells] <- NA
  
  # Generate a set of initiation probability rasters ---------------------------
  
  initiationProbRasterList <- list()
  
  for (i in seq_len(iterations)) {
    ## Generate non-initiation regions -----------------------------------------
    
    # Sample points in areas that fit initiation ranges but recorded no landslides
    noninitiationPoints <- terra::spatSample(
      initiationRaster, size = noninitiationPointsCount, na.rm = TRUE, as.points = TRUE
    )
    
    # Create buffer regions around non-initiation points
    noninitiationPolys <- terra::buffer(noninitiationPoints, width = bufferRadius)
    
    ## Create training dataset -------------------------------------------------
    
    # Extract variable values within non-initiation regions
    noninitiationValues <- terra::extract(varsRaster, noninitiationPolys)
    noninitiationValues <- noninitiationValues[-1] # Remove "ID" column
    
    # Assign a classification value to the initiation and non-initiation entries
    initiationValues$class <- rep("initiation", nrow(initiationValues))
    noninitiationValues$class <- rep("non-initiation", nrow(noninitiationValues))
    
    # Combine initiation and non-initiation entries into a single dataset
    trainingData <- rbind(initiationValues, noninitiationValues)
    
    # Filter out entries with NA values
    trainingData <- na.omit(trainingData)
    
    # Factor the classification variable values
    trainingData$class <- factor(trainingData$class)
    
    ## Create random forest model ----------------------------------------------
    
    rfModel <- randomForest::randomForest(
      formula = class ~ .,
      data = trainingData
    )
    
    ## Generate probability raster ---------------------------------------------
    
    initiationProbRaster <- terra::predict(
      varsRaster,
      rfModel,
      na.rm = TRUE,
      type = "prob"
    )[["initiation"]]
    
    initiationProbRasterList[[i]] <- initiationProbRaster
  }
  
  # Summarize initiation probability rasters -----------------------------------
  
  if (length(initiationProbRasterList) > 1) {
    # Combine probability rasters into a single multi-layer raster
    initiationProbRaster <- terra::rast(initiationProbRasterList)
    
    # Generate summary initiation probability rasters
    avgInitionProbRaster <- terra::app(initiationProbRaster, fun = "mean")
    minInitionProbRaster <- terra::app(initiationProbRaster, fun = "min")
    maxInitionProbRaster <- terra::app(initiationProbRaster, fun = "max")
    
    # Save summary rasters
    terra::writeRaster(avgInitionProbRaster, paste0(outputDir, "/prob_avg.tif"), overwrite = TRUE)
    terra::writeRaster(minInitionProbRaster, paste0(outputDir, "/prob_min.tif"), overwrite = TRUE)
    terra::writeRaster(maxInitionProbRaster, paste0(outputDir, "/prob_max.tif"), overwrite = TRUE)
  } else {
    # Save single probability raster
    terra::writeRaster(initiationProbRasterList[[1]], paste0(outputDir, "/prob.tif"), overwrite = TRUE)
  }
  
  logMsg("Wrote summary initiation probability rasters")
  
  # Return ---------------------------------------------------------------------
  
  return(out_params)
  
}

if (FALSE) {
  
  tool_exec(
    in_params = list(
      referenceRasterFile = "C:/Work/netmapdata/Scottsburg/elev_scottsburg.flt",
      varRasterFiles <- list(
        "C:/Work/netmapdata/Scottsburg/grad_30.tif",
        "C:/Work/netmapdata/Scottsburg/plan_15.tif"
      ),
      initiationPointsFile = "C:/Work/netmapdata/Scottsburg/Scottsburg_Upslope.shp",
      bufferRadius = 20,
      initiationLimitScaler = 1.05,
      noninitiationPointsCount = 50,
      iterations = 2,
      outputDir = "C:/Work/netmapdata/Scottsburg"
    ),
    out_params = list()
  )
  
  # Test in Scottsburg (WORK2)
  tool_exec(
    in_params = list(
      referenceRasterFile = "E:/NetmapData/Scottsburg/elev_scottsburg.flt",
      varRasterFiles <- list(
        "E:/NetmapData/Scottsburg/grad_30.tif",
        "E:/NetmapData/Scottsburg/plan_30.tif"
      ),
      initiationPointsFile = "E:/NetmapData/Scottsburg/Scottsburg_Upslope.shp",
      bufferRadius = 20,
      initiationLimitMultiplier = 1.05,
      noninitiationPointsCount = 50,
      iterations = 2,
      outputDir = "E:/NetmapData/Scottsburg"
    ),
    out_params = list()
  )
  
}
