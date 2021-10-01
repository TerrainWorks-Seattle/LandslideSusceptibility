tool_exec <- function(in_params, out_params) {
  
  # Set parameters -------------------------------------------------------------
  
  # Input
  referenceRasterFile <- in_params[[1]]
  varRasterFiles <- in_params[[2]]
  initiationPointsFile <- in_params[[3]]
  bufferRadius <- in_params[[4]]
  initiationLimitPercent <- in_params[[5]]
  k <- in_params[[6]]
  generateProbabilityRasters <- in_params[[7]]
  outputDir <- in_params[[8]]
  
  # Set up logging -------------------------------------------------------------
  
  logFilename <- "landslide_susceptibility_validation_log.txt"
  file.create(paste0(outputDir, "/", logFilename))
  
  logObj <- function(obj) {
    capture.output(obj, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  logMsg <- function(msg) {
    cat(msg, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
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
  
  # Calculate variables' initiation range --------------------------------------
  
  # Extract all variable values from initiation buffers
  initiationValues <- terra::extract(varsRaster, initiationBuffers)
  
  # For each initiation buffer, find each variable's maximum value
  regionMaxVarValues <- aggregate(. ~ ID, data = initiationValues, max)
  regionMaxVarValues <- regionMaxVarValues[-1] # Remove "ID" column
  initiationValues <- initiationValues[-1] # Remove the "ID" column
  
  # Find the min and max of the maximum variable values in each buffer 
  initiationMinValues <- lapply(regionMaxVarValues, min)
  initiationMaxValues <- lapply(regionMaxVarValues, max)
  
  # Slightly expand initiation range
  initiationLimitRatio <- initiationLimitPercent / 100  
  
  initiationMinValues <- lapply(initiationMinValues, function(x) x - initiationLimitRatio * x)
  initiationMaxValues <- lapply(initiationMaxValues, function(x) x + initiationLimitRatio * x)
  
  # Define matrix to hold each variable's initiation value range
  initiationRange <- matrix(rep(NA, 4), nrow = 2)
  colnames(initiationRange) <- c("min", "max")
  rownames(initiationRange) <- names(varsRaster)
  
  # Populate matrix with range limits
  for (varName in names(varsRaster)) {
    initiationRange[varName, "min"] <- initiationMinValues[[varName]]
    initiationRange[varName, "max"] <- initiationMaxValues[[varName]]
  }
  
  # Create initiation mask -----------------------------------------------------
  
  # Define a raster which will NA-out any cells with variable values outside of 
  # their initiation ranges
  initiationRaster <- terra::copy(varsRaster)
  
  for (varName in names(initiationRaster)) {
    varRaster <- initiationRaster[[varName]]
    
    # Get variable value limits
    minInitiationValue <- initiationMinValues[[varName]]
    maxInitiationValue <- initiationMaxValues[[varName]]
    
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
  initiationCells <- terra::extract(initiationRaster, initiationBuffers, cells = TRUE)[["cell"]]
  initiationRaster[initiationCells] <- NA
  
  # Generate non-initiation regions --------------------------------------------
  
  # NOTE: terra::spatSample() sometimes generates less than the requested 
  # number of points if the sample raster has a lot of NAs. This is remedied 
  # by iteratively requesting a larger and larger number of points until
  # enough have been generated, then subsetting those.
  
  hasGeneratedEnough <- FALSE
  pointsPerClass <- length(initiationPoints)
  requestCount <- pointsPerClass
  
  while (!hasGeneratedEnough) {
    # Sample points in areas that fit initiation ranges but recorded no landslides
    noninitiationPoints <- terra::spatSample(
      initiationRaster,
      size = requestCount,
      na.rm = TRUE,
      as.points = TRUE,
      warn = FALSE
    )
    
    # Exit once enough non-initiation points have been generated
    if (length(noninitiationPoints) >= pointsPerClass) {
      noninitiationPoints <- noninitiationPoints[seq_along(initiationPoints)]
      hasGeneratedEnough <- TRUE
    }
    
    requestCount <- requestCount * 10
  }
  
  # Create a buffer around each non-initiation point
  noninitiationBuffers <- if (bufferRadius > 0) {
    terra::buffer(noninitiationPoints, width = bufferRadius)
  } else {
    noninitiationPoints
  }
  
  # Create buffer subsets ------------------------------------------------------
  
  totalBuffers <- length(initiationBuffers) + length(noninitiationBuffers)
  testingBuffersPerSubset <- floor(totalBuffers / k)
  
  testingFreeInitiationIndices <- seq_along(initiationBuffers)
  testingFreeNoninitiationIndices <- seq_along(noninitiationBuffers)
  
  testingSets <- list()
  for (i in 1:k) {
    
    # How many of each buffer type should be in the testing set?
    # NOTE: An equal number for now, otherwise there could potentially be zero
    # of one type for an iteration
    testingInitiationCount <- length(initiationBuffers) / k / 2
    testingNoninitiationCount <- testingBuffersPerSubset - testingInitiationCount
    
    # Choose the initiation buffers for this testing set
    if (length(testingFreeInitiationIndices) > 0) {
      testingInitiationIndices <- sample(testingFreeInitiationIndices, size = min(7, length(testingFreeInitiationIndices)))
      testingFreeInitiationIndices <- testingFreeInitiationIndices[! testingFreeInitiationIndices %in% testingInitiationIndices]
    } else {
      testingInitiationIndices <- NA
    }
    
    # Choose the non-initiation buffers for this testing set
    if (length(testingFreeNoninitiationIndices) > 0) {
      testingNoninitiationIndices <- sample(testingFreeNoninitiationIndices, size = min(7, length(testingFreeNoninitiationIndices)))
      testingFreeNoninitiationIndices <- testingFreeNoninitiationIndices[! testingFreeNoninitiationIndices %in% testingNoninitiationIndices]
    } else {
      testingNoninitiationIndices <- NA
    }
    
    # Store the test set
    testingSets[[i]] <- list(
      initiation = testingInitiationIndices,
      noninitiation = testingNoninitiationIndices
    )
    
  }
  
  # Perform cross-validation ---------------------------------------------------
  
  # Place to store iteration auc values
  iterationsAucValues <- c()
  
  # Place to store iteration error rates
  iterationsErrorRates <- data.frame(
    rep(NA, k),
    rep(NA, k),
    rep(NA, k)
  )
  names(iterationsErrorRates) <- c("OOB", "initiation", "non-initiation")
  
  for (i in seq_along(testingSets)) {
    logMsg(paste0("Model ", formatC(i, width = 2, format = "d", flag = "0"), 
                  " -------------------------------------------\n"))
    
    ## Create training dataset -------------------------------------------------
    
    testingSet <- testingSets[[i]]
    
    # Get training buffers (all buffers not in testing set)
    trainingInitiationIndices <- setdiff(seq_along(initiationBuffers), testingSet$initiation)
    trainingNoninitiationIndices <- setdiff(seq_along(noninitiationBuffers), testingSet$noninitiation)
    
    trainingInitiationBuffers <- initiationBuffers[trainingInitiationIndices]
    trainingNoninitiationBuffers <- noninitiationBuffers[trainingNoninitiationIndices]
    
    # Extract values from training initiation buffers
    trainingInitiationValues <- terra::extract(varsRaster, trainingInitiationBuffers)
    trainingNoninitiationValues <- terra::extract(varsRaster, trainingNoninitiationBuffers)
    
    # Assign a classification value to each entry
    trainingInitiationValues$class <- rep("initiation", nrow(trainingInitiationValues))
    trainingNoninitiationValues$class <- rep("non-initiation", nrow(trainingNoninitiationValues))
    
    # Combine initiation and non-initiation entries into a single dataset
    trainingData <- rbind(trainingInitiationValues, trainingNoninitiationValues)
    
    # Filter out entries with NA values
    trainingData <- na.omit(trainingData)
    
    # Factor the classification variable values
    trainingData$class <- factor(trainingData$class)
    
    # Remove "ID" column
    trainingData$ID <- NULL
    
    ## Create model ------------------------------------------------------------
    
    # Create random forest model
    rfModel <- randomForest::randomForest(
      formula = class ~ .,
      data = trainingData
    )
    
    # Log model error rates
    logMsg("Error rates:\n")
    logObj(rfModel$err.rate[rfModel$ntree,])
    logMsg("\n")
    
    # Log model confusion matrix
    logMsg("Model confusion matrix:\n")
    logObj(rfModel$confusion)
    logMsg("\n")
    
    ## Generate probability raster --------------------------------------------
    
    if (generateProbabilityRasters) {
      # Generate probability raster
      initiationProbRaster <- terra::predict(
        varsRaster,
        rfModel,
        na.rm = TRUE,
        type = "prob"
      )[["initiation"]]
      
      # Save raster
      terra::writeRaster(initiationProbRaster, paste0(outputDir, "/prob_", i, ".tif"))
    }
    
    ## Create testing dataset ----------------------------------------------
    
    # Get testing buffers
    testingInitiationIndices <- testingSet$initiation
    testingNoninitiationIndices <- testingSet$noninitiation
    
    testingInitiationBuffers <- initiationBuffers[testingInitiationIndices]
    testingNoninitiationBuffers <- noninitiationBuffers[testingNoninitiationIndices]
    
    # Extract values from testing initiation buffers
    testingInitiationValues <- terra::extract(varsRaster, testingInitiationBuffers)
    testingNoninitiationValues <- terra::extract(varsRaster, testingNoninitiationBuffers)
    
    # Assign a classification value to each entry
    testingInitiationValues$class <- rep("initiation", nrow(testingInitiationValues))
    testingNoninitiationValues$class <- rep("non-initiation", nrow(testingNoninitiationValues))
    
    # Combine initiation and non-initiation entries into a single dataset
    testingData <- rbind(testingInitiationValues, testingNoninitiationValues)
    
    # Filter out entries with NA values
    testingData <- na.omit(testingData)
    
    # Factor the classification variable values
    testingData$class <- factor(testingData$class)
    
    # Remove "ID" column
    testingData$ID <- NULL
    
    ## Run model on testing data ----------------------------------------------
    
    # Predict the class type of each test dataset entry
    predictedClass <- predict(
      rfModel,
      type = "response",
      newdata = testingData
    )
    
    # Log test dataset confusion matrix
    logMsg("Testing confusion matrix:")
    logObj(table(predictedClass, testingData$class))
    logMsg("\n")
    
    ## Run model on testing data -----------------------------------------------
    
    # Calculate initiation probability for each test entry
    initiationProb <- predict(
      rfModel,
      type = "prob",
      newdata = testingData
    )[,"initiation"]
    
    # Calculate ROC stats
    rocStats <- TerrainWorksUtils::calcRocStats(
      classes = testingData$class,
      probs = initiationProb,
      "initiation",
      "non-initiation"
    )
    
    # Log AUC value
    auc <- rocStats$auc@y.values[[1]]
    logMsg(paste0("AUC: ", auc, "\n\n"))
    
    ## Record iteration statistics ---------------------------------------------
    
    iterationsAucValues <- c(iterationsAucValues, auc)
    iterationsErrorRates[i,] <- as.data.frame(rfModel$err.rate)[rfModel$ntree,]
  }
  
  # Summarize iterations -------------------------------------------------------
  
  logMsg("Summary ------------------------------------------\n")
  
  # Calculate AUC summary stats
  aucMin <- min(iterationsAucValues)
  aucMax <- max(iterationsAucValues)
  aucRange <- diff(c(aucMin, aucMax))
  aucStdev <- sd(iterationsAucValues)
  
  aucMatrix <- matrix(c(aucMin, aucMax, aucRange, aucStdev), ncol = 1)
  rownames(aucMatrix) <- c("Min", "Max", "Range", "Stdev")
  
  # Log standard deviation of AUC values
  logMsg("AUC stats:\n")
  logObj(aucMatrix)
  logMsg("\n")
  
  # Calculate error rates stats
  errorRateMin <- as.data.frame(lapply(iterationsErrorRates, min, na.rm = TRUE))
  errorRateMax <- as.data.frame(lapply(iterationsErrorRates, max, na.rm = TRUE))
  errorRateRange <- as.data.frame(lapply(iterationsErrorRates, function(x) {
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  }))
  errorRateStdev <- as.data.frame(lapply(iterationsErrorRates, sd, na.rm = TRUE))
  
  errorRateDf <- rbind(
    errorRateMin,
    errorRateMax,
    errorRateRange,
    errorRateStdev
  )
  
  errorRateMatrix <- as.matrix(errorRateDf)
  rownames(errorRateMatrix) <- c("Min", "Max", "Range", "Stdev")
  
  # Log error rates stats
  logMsg("Error rate stats:\n")
  logObj(errorRateMatrix)
  
  # Return ---------------------------------------------------------------------
  
  return(out_params)
  
}


if (FALSE) {
  
  # Run in Scottsburg (BIGLAPTOP)
  tool_exec(
    in_params = list(
      referenceRasterFile = "C:/Work/netmapdata/Scottsburg/elev_scottsburg.flt",
      varRasterFiles = list(
        "C:/Work/netmapdata/Scottsburg/grad_30.tif",
        "C:/Work/netmapdata/Scottsburg/plan_30.tif"
      ),
      initiationPointsFile = "C:/Work/netmapdata/Scottsburg/Scottsburg_Upslope.shp",
      bufferRadius = 20,
      initiationLimitPercent = 20,
      k = 5,
      generateProbabilityRasters = FALSE,
      outputDir = "C:/Work/netmapdata/Scottsburg"
    ),
    out_params = list()
  )
  
}
