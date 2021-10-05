# This script provides a way to assess how much the choice of initiation points
# affects the accuracy of a landslide susceptibility prediction model.
# k iterations of the model are created, each trained with a unique set of
# initiation points and a static set of non-initiation points. A buffer of given
# size is drawn around each point and all cell values within are used as model 
# inputs.

tool_exec <- function(in_params, out_params) {
  
  # Set parameters -------------------------------------------------------------
  
  # Input
  referenceRasterFile <- in_params[[1]]
  varRasterFiles <- in_params[[2]]
  initiationPointsFile <- in_params[[3]]
  bufferRadius <- in_params[[4]]
  initiationLimitPercent <- in_params[[5]]
  noninitiationRatio <- in_params[[6]]
  k <- in_params[[7]]
  generateProbabilityRasters <- in_params[[8]]
  outputDir <- in_params[[9]]
  
  # Set up logging -------------------------------------------------------------
  
  logFilename <- "initiation_sensitivity_full_buffer.txt"
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
  
  # Generate non-initiation buffers --------------------------------------------
  
  # NOTE: terra::spatSample() sometimes generates less than the requested 
  # number of points if the sample raster has a lot of NAs. This is remedied 
  # by repeatedly requesting a larger and larger number of points until
  # enough have been generated, then subsetting those.
  
  desiredNoninitiationCount <- ceiling(length(initiationPoints) * noninitiationRatio)
  hasGeneratedEnough <- FALSE
  requestCount <- desiredNoninitiationCount
  
  while (!hasGeneratedEnough) {
    # Sample points anywhere that fits initiation conditions but recorded no landslides
    noninitiationPoints <- terra::spatSample(
      initiationRaster,
      size = requestCount,
      na.rm = TRUE,
      as.points = TRUE,
      warn = FALSE
    )
    
    # Exit once enough non-initiation points have been generated
    if (length(noninitiationPoints) >= desiredNoninitiationCount) {
      noninitiationPoints <- noninitiationPoints[seq_len(desiredNoninitiationCount)]
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
  
  # Create initiation buffer sets ----------------------------------------------
  
  # Calculate how many initiation buffers there should be for testing each iteration
  testingInitiationBuffersPerIteration <- floor(length(initiationBuffers) / k)
  
  # Create k unique sets of initiation buffers to use for training and testing
  testingSets <- list()
  testingFreeInitiationIndices <- seq_along(initiationBuffers)
  for (i in 1:k) {
    testingInitiationIndices <- sample(testingFreeInitiationIndices, size = testingInitiationBuffersPerIteration)
    testingFreeInitiationIndices <- testingFreeInitiationIndices[!(testingFreeInitiationIndices %in% testingInitiationIndices)]
    testingSets[[i]] <- testingInitiationIndices
  }
  
  # Create non-initiation buffer sets ------------------------------------------
  
  # Create static sets of training and testing non-initiation buffers to use every iteration
  testingNoninitiationCount <- floor(testingInitiationBuffersPerIteration * noninitiationRatio)
  testingNoninitiationIndices <- sample(seq_along(noninitiationBuffers), size = testingNoninitiationCount)
  trainingNoninitiationIndices <- setdiff(seq_along(noninitiationBuffers), testingNoninitiationIndices)
  
  testingNoninitiationBuffers <- noninitiationBuffers[testingNoninitiationIndices]
  trainingNoninitiationBuffers <- noninitiationBuffers[trainingNoninitiationIndices]
  
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
                  " -------------------------------------------\n\n"))
    
    ## Create training dataset -------------------------------------------------
    
    # Get training initiation buffers
    trainingInitiationIndices <- setdiff(seq_along(initiationBuffers), testingSets[[i]])
    trainingInitiationBuffers <- initiationBuffers[trainingInitiationIndices]
    
    # Extract values from training buffers
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
    logMsg("ERROR RATES:\n")
    errorRateDf <- data.frame(rfModel$err.rate[rfModel$ntree,])
    colnames(errorRateDf) <- "error rate"
    logObj(errorRateDf)
    logMsg("\n")
    
    # Log model confusion matrix
    logMsg("MODEL CONFUSION MATRIX:\n")
    logObj(rfModel$confusion)
    logMsg("\n")
    
    ## Generate probability raster ---------------------------------------------
    
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
    
    ## Create testing dataset --------------------------------------------------
    
    # Get testing buffers
    testingInitiationIndices <- testingSets[[i]]
    testingInitiationBuffers <- initiationBuffers[testingInitiationIndices]
    
    # Extract values from testing buffers
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
    
    ## Run model on testing data -----------------------------------------------
    
    # Predict the class type of each test dataset entry
    predictedClass <- predict(
      rfModel,
      type = "response",
      newdata = testingData
    )
    
    # Log test dataset confusion matrix
    logMsg("TESTING CONFUSION MATRIX:")
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
    logMsg(paste0("AUC: ", round(auc, digits = 7), "\n\n"))
    
    ## Record iteration statistics ---------------------------------------------
    
    iterationsAucValues <- c(iterationsAucValues, auc)
    iterationsErrorRates[i,] <- as.data.frame(rfModel$err.rate)[rfModel$ntree,]
  }
  
  # Summarize iterations -------------------------------------------------------
  
  logMsg("SUMMARY --------------------------------------------\n\n")
  
  # Calculate AUC summary stats
  aucMin <- min(iterationsAucValues)
  aucMax <- max(iterationsAucValues)
  aucRange <- diff(c(aucMin, aucMax))
  aucStdev <- sd(iterationsAucValues)
  
  aucMatrix <- matrix(c(aucMin, aucMax, aucRange, aucStdev), ncol = 1)
  aucMatrix <- t(aucMatrix)
  colnames(aucMatrix) <- c("Min", "Max", "Range", "Stdev")
  rownames(aucMatrix) <- "AUC"
  
  # Log standard deviation of AUC values
  logMsg("AUC:\n")
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
  
  errorRateDf <- t(errorRateDf)
  
  errorRateMatrix <- as.matrix(errorRateDf)
  colnames(errorRateMatrix) <- c("min", "max", "range", "stdev")
  
  # Log error rates stats
  logMsg("ERROR RATES:\n")
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
      initiationLimitPercent = 10,
      noninitiationRatio = 1.5,
      k = 5,
      generateProbabilityRasters = FALSE,
      outputDir = "C:/Work/netmapdata/Scottsburg"
    ),
    out_params = list()
  )
  
 }
 