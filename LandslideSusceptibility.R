tool_exec <- function(in_params, out_params) {
  # Helper functions -----------------------------------------------------------
  
  # Creates a matrix that holds the min/max initiation limits of each 
  # explanatory variable.
  createInitiationRange <- function(
    varsRaster,              # A SpatRaster of explanatory variables
    initiationBuffers,       # A SpatVector of initiation buffers
    initiationRangeExpansion # Percent to reduce each range min and increase each range max by 
  ) {
    # Extract all variable values from initiation buffers
    initiationValues <- terra::extract(varsRaster, initiationBuffers)
    
    # For each initiation buffer, find each variable's maximum value
    regionMaxVarValues <- aggregate(. ~ ID, data = initiationValues, max)
    regionMaxVarValues$ID <- NULL # Remove "ID" column
    initiationValues$ID <- NULL # Remove the "ID" column
    
    # Find the min and max of the maximum variable values in each buffer 
    initiationMinValues <- lapply(regionMaxVarValues, min)
    initiationMaxValues <- lapply(regionMaxVarValues, max)
    
    # Slightly expand initiation range
    r <- initiationRangeExpansion / 100  
    initiationMinValues <- lapply(initiationMinValues, function(x) x - r * x)
    initiationMaxValues <- lapply(initiationMaxValues, function(x) x + r * x)
    
    # Create a matrix that holds each variable's initiation range
    initiationRange <- matrix(rep(NA, 4), nrow = 2)
    colnames(initiationRange) <- c("min", "max")
    rownames(initiationRange) <- names(varsRaster)
    for (varName in names(varsRaster)) {
      initiationRange[varName, "min"] <- initiationMinValues[[varName]]
      initiationRange[varName, "max"] <- initiationMaxValues[[varName]]
    }
    
    return(initiationRange)
  }
  
  # Creates a raster with NA cells anywhere the input variable rasters don't
  # fall within their respective initiation ranges or are inside initiation
  # sites.
  createNoninitiationRaster <- function(
    varsRaster,       # A SpatRaster of explanatory variables
    initiationRange,  # A matrix that holds the min/max initiation limits of each explanatory variable 
    initiationBuffers # A SpatVector of initiation buffers
  ) {
    initiationRaster <- terra::copy(varsRaster)
    for (varName in names(initiationRaster)) {
      varRaster <- initiationRaster[[varName]]
      
      # Get variable value limits
      minInitiationValue <- initiationRange[varName, "min"]
      maxInitiationValue <- initiationRange[varName, "max"]
      
      # NA-out cells with values outside variable initiation range
      varInitiationRaster <- terra::app(varRaster, function(x) {
        ifelse(x < minInitiationValue | x > maxInitiationValue, NA, x)
      })
      
      # Update the raster in the input raster stack
      initiationRaster[[varName]] <- varInitiationRaster
    }
    
    # NA-out cells with variable values outside their initiation range
    initiationRaster <- terra::app(initiationRaster, fun = "prod")
    
    # NA-out cells within initiation buffers
    initiationCells <- terra::extract(initiationRaster, initiationBuffers, cells = TRUE)$cell
    initiationRaster[initiationCells] <- NA
    
    return(initiationRaster)
  }
  
  # Generates a SpatVector of non-initiation buffers.
  generateNoninitiationBuffers <- function(
    noninitiationBuffersCount, # The number of non-initiation buffers to generate
    noninitiationRegion,       # A SpatRaster delineating where non-initiation buffers can be placed
    bufferRadius               # The radius of each non-initiation buffer
  ) {
    # NOTE: terra::spatSample() sometimes generates less than the requested 
    # number of points if the sample raster has a lot of NAs. This is remedied 
    # by repeatedly requesting a larger and larger number of points until
    # enough have been generated, then subsetting those for the correct amount.
    
    currentRequest <- noninitiationBuffersCount
    hasGeneratedEnough <- FALSE
    
    while (!hasGeneratedEnough) {
      # Sample points anywhere that fits initiation conditions but recorded no landslides
      noninitiationPoints <- terra::spatSample(
        noninitiationRegion,
        size = currentRequest,
        na.rm = TRUE,
        as.points = TRUE,
        warn = FALSE
      )
      
      # Test if enough non-initiation points have been generated
      if (length(noninitiationPoints) >= noninitiationBuffersCount) {
        # If so, subset and exit loop
        noninitiationPoints <- noninitiationPoints[seq_len(noninitiationBuffersCount)]
        hasGeneratedEnough <- TRUE
      } else {
        # If not, double the next request
        currentRequest <- currentRequest * 2
      }
    }
    
    # Create a buffer around each non-initiation point
    noninitiationBuffers <- if (bufferRadius > 0) {
      terra::buffer(noninitiationPoints, width = bufferRadius)
    } else {
      noninitiationPoints
    }
    
    return(noninitiationBuffers)
  }
  
  # Create a dataset containing variable values and classes of cells within 
  # initiation and non-initiation buffers.
  extractBufferValues <- function(
    varsRaster,            # A SpatRaster of explanatory variables
    initiationBuffers,     # A SpatVector of initiation buffers
    noninitiationBuffers,  # A SpatVector of non-initiation buffers,
    bufferExtractionMethod # Method to use for extracting buffer values 
  ) {
    # By default, extract all values from initiation and non-initiation buffers
    initiationValues <- terra::extract(varsRaster, initiationBuffers)
    noninitiationValues <- terra::extract(varsRaster, noninitiationBuffers)
    
    # Subset if a different buffer extraction method was requested
    if (bufferExtractionMethod == "center cell") {
      # Find center points of buffers
      initiationCenters <- terra::centroids(initiationBuffers)
      noninitiationCenters <- terra::centroids(noninitiationBuffers)
      
      # Get center point coordinates
      initiationCoords <- terra::geom(initiationCenters)[,c("x","y")]
      noninitiationCoords <- terra::geom(noninitiationCenters)[,c("x","y")]
      
      # Extract values from cells containing center points
      initiationValues <- terra::extract(varsRaster, initiationCoords)
      noninitiationValues <- terra::extract(varsRaster, noninitiationCoords)
    } else if (bufferExtractionMethod == "max gradient cell") {
      initiationValues <- aggregateBufferValues(initiationValues, "grad", max)
      noninitiationValues <- aggregateBufferValues(noninitiationValues, "grad", max)
    } else if (bufferExtractionMethod == "max plan cell") {
      initiationValues <- aggregateBufferValues(initiationValues, "plan", max)
      noninitiationValues <- aggregateBufferValues(noninitiationValues, "plan", max)
    }
    
    # Assign a classification value to each entry
    initiationValues$class <- rep("initiation", nrow(initiationValues))
    noninitiationValues$class <- rep("non-initiation", nrow(noninitiationValues))
    
    # Combine initiation and non-initiation entries into a single dataset
    dataset <- rbind(initiationValues, noninitiationValues)
    
    # Remove the "ID" column
    dataset$ID <- NULL
    
    # Factor the classification variable values
    dataset$class <- factor(dataset$class)
    
    # Filter out entries with NA values
    dataset <- na.omit(dataset)
    
    return(dataset)
  }
  
  # Groups values by buffer and, for each buffer, keeps the entry with the "fun"
  # (min/max) variable value.
  aggregateBufferValues <- function(
    bufferValues, # Values extracted from buffers, including an "ID" column
    varName,      # Name pattern of the variable to aggregate by
    fun           # Function to aggregate buffer values by
  ) {
    # Determine which variable to aggregate by based on provided name pattern
    gradientVarName <- names(varsRaster)[grepl("grad", names(varsRaster))][1]
    
    # Formula to group entries by buffer and return requested variable value
    fm <- as.formula(paste(gradientVarName, "~", "ID"))
    
    # For each buffer, keep the entry with the "fun" variable value
    aggregatedValues <- merge(
      aggregate(fm, max, data = bufferValues),
      bufferValues
    )
    
    return(aggregatedValues)
  }
  
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
