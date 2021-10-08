tool_exec <- function(in_params, out_params) {
  
  # Load helper functions ------------------------------------------------------
  
  source("./shared/createInitiationRange.R")
  source("./shared/createNoninitiationRaster.R")
  source("./shared/generateNoninitiationBuffers.R")
  source("./shared/extractBufferValues.R")
  
  # Set parameters -------------------------------------------------------------
  
  # Input
  referenceRasterFile      <- in_params[[1]]
  varRasterFiles           <- in_params[[2]]
  initiationPointsFile     <- in_params[[3]]
  noninitiationRatio       <- in_params[[4]]
  bufferRadius             <- in_params[[5]]
  bufferExtractionMethod   <- in_params[[6]]
  initiationRangeExpansion <- in_params[[7]]
  iterations               <- in_params[[8]]
  outputDir                <- in_params[[9]]
  
  # Set up logging -------------------------------------------------------------
  
  logFilename <- "landslide_susceptibility.txt"
  file.create(paste0(outputDir, "/", logFilename))
  
  logObj <- function(obj) {
    capture.output(obj, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  logMsg <- function(msg) {
    cat(msg, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  logMsg("INPUT PARAMETERS:\n")
  logMsg(paste0("  Reference raster: ", referenceRasterFile, "\n"))
  logMsg("  Explanatory variable rasters:\n")
  for (i in seq_along(varRasterFiles)) { logMsg(paste0("    [", i, "] ", varRasterFiles[[i]], "\n")) }
  logMsg(paste0("  Initiation points: ", initiationPointsFile, "\n"))
  logMsg(paste0("  Non-initiation points ratio: ", noninitiationRatio, "\n"))
  logMsg(paste0("  Buffer radius: ", bufferRadius, "\n"))
  logMsg(paste0("  Buffer extraction method: ", bufferExtractionMethod, "\n"))
  logMsg(paste0("  Initiation range expansion: ", initiationRangeExpansion, "%\n"))
  logMsg(paste0("  Iterations: ", iterations, "\n"))
  logMsg(paste0("  Output directory: ", outputDir, "\n"))
  logMsg("\n")
  
  # Load rasters ---------------------------------------------------------------
  
  # Load reference raster
  referenceRaster <- terra::rast(referenceRasterFile)
  
  # Load explanatory variable rasters
  varRasterList <- lapply(varRasterFiles, function(file) terra::rast(file))
  
  # Align variable rasters
  varRasterList <- TerrainWorksUtils::alignRasters(referenceRaster, varRasterList)
  
  # Combine variable rasters into a single multi-layer raster
  varsRaster <- terra::rast(varRasterList)
  
  # Create initiation buffers --------------------------------------------------
  
  # Load initiation points
  initiationPoints <- terra::vect(initiationPointsFile)
  initiationPoints <- terra::project(initiationPoints, referenceRaster)
  
  # Create a buffer around each initiation point
  initiationBuffers <- if (bufferRadius > 0) {
    terra::buffer(initiationPoints, width = bufferRadius)
  } else {
    initiationPoints
  }
  
  # Generate non-initiation buffers --------------------------------------------
  
  # Calculate initiation range for each variable
  initiationRange <- createInitiationRange(
    varsRaster,
    initiationBuffers,
    initiationRangeExpansion
  )
  
  # Log initiation range matrix
  logMsg("INITIATION RANGES:\n")
  logObj(initiationRange)
  logMsg("\n")
  
  # The region where non-initiation buffers can be sampled from: regions that 
  # meet initiation conditions but recorded no landslides
  noninitiationRaster <- createNoninitiationRaster(
    varsRaster,
    initiationRange,
    initiationBuffers
  )
  
  # Generate a set of initiation probability rasters ---------------------------
  
  # Place to store initation probability rasters
  initiationProbRasterList <- list()
  
  for (i in seq_len(iterations)) {
    ## Generate non-initiation buffers -----------------------------------------
    
    # Determine how many to generate
    noninitiationBuffersCount <- ceiling(length(initiationPoints) * noninitiationRatio)
    
    # Generate non-initiation buffers
    noninitiationBuffers <- generateNoninitiationBuffers(
      noninitiationBuffersCount,
      noninitiationRaster,
      bufferRadius
    )
      
    ## Create model ------------------------------------------------------------
    
    # Create training dataset
    trainingData <- extractBufferValues(
      varsRaster,
      initiationBuffers,
      noninitiationBuffers,
      bufferExtractionMethod
    )
    
    # Train a new random forest model
    rfModel <- randomForest::randomForest(
      formula = class ~ .,
      data = trainingData
    )
    
    logMsg(paste0("Model ", formatC(i, width = 2, format = "d", flag = "0"), 
                  " ------------------------------------------------\n\n"))
    
    logMsg("MODEL ERROR RATES:\n")
    errorRateDf <- data.frame(rfModel$err.rate[rfModel$ntree,])
    colnames(errorRateDf) <- "error rate"
    logObj(errorRateDf)
    logMsg("\n")
    
    # Log model confusion matrix
    logMsg("MODEL CONFUSION MATRIX:\n")
    rownames(rfModel$confusion) <- c("true initiation", "true non-initiation")
    logObj(rfModel$confusion)
    logMsg("\n")
    
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
  
  # Return ---------------------------------------------------------------------
  
  return(out_params)
  
}


if (FALSE) {
  
  tool_exec(
    in_params = list(
      referenceRasterFile = "C:/Work/netmapdata/Scottsburg/elev_scottsburg.flt",
      varRasterFiles <- list(
        "C:/Work/netmapdata/Scottsburg/grad_30.tif",
        "C:/Work/netmapdata/Scottsburg/plan_30.tif"
      ),
      initiationPointsFile = "C:/Work/netmapdata/Scottsburg/Scottsburg_Upslope.shp",
      noninitiationRatio = 1.5,
      bufferRadius = 20,
      bufferExtractionMethod = "center cell",
      initiationRangeExpansion = 10,
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
      noninitiationRatio = 1.5,
      bufferRadius = 20,
      bufferExtractionMethod = "center cell",
      initiationRangeExpansion = 10,
      iterations = 2,
      outputDir = "E:/NetmapData/Scottsburg"
    ),
    out_params = list()
  )
  
}
