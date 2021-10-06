# This script provides a way to assess how much the choice of initiation points
# affects the accuracy of a landslide susceptibility prediction model.
# k iterations of the model are created, each trained with a unique set of
# initiation points and a static set of non-initiation points. A buffer of given
# size is drawn around each point and all cell values within are used as model 
# inputs.

tool_exec <- function(in_params, out_params) {
  
  # Helper functions -----------------------------------------------------------
  
  # Creates a matrix that holds the min/max initiation limits of each 
  # explanatory variable.
  createInitiationRange <- function(
    varsRaster,            # A SpatRaster of explanatory variables
    initiationBuffers,     # A SpatVector of initiation buffers
    initiationLimitPercent # Percent to reduce each range min and increase each range max by 
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
    initiationLimitRatio <- initiationLimitPercent / 100  
    initiationMinValues <- lapply(initiationMinValues, function(x) x - initiationLimitRatio * x)
    initiationMaxValues <- lapply(initiationMaxValues, function(x) x + initiationLimitRatio * x)
    
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
    if (bufferExtractionMethod == "steepest cell") {
      # Group by buffer
      # Keep entry with the max grad value of each group
    } else if (bufferExtractionMethod == "most convergent cell") {
      # Group by buffer
      # Keep entry with the max plan value of each group
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
  
  # Set parameters -------------------------------------------------------------
  
  # Input
  referenceRasterFile <- in_params[[1]]
  varRasterFiles <- in_params[[2]]
  initiationPointsFile <- in_params[[3]]
  noninitiationRatio <- in_params[[4]]
  bufferRadius <- in_params[[5]]
  bufferExtractionMethod <- in_params[[6]]
  initiationLimitPercent <- in_params[[7]]
  k <- in_params[[8]]
  generateProbabilityRasters <- in_params[[9]]
  outputDir <- in_params[[10]]
  
  # Set up logging -------------------------------------------------------------
  
  logFilename <- "initiation_sensitivity.txt"
  file.create(paste0(outputDir, "/", logFilename))
  
  logObj <- function(obj) {
    capture.output(obj, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  logMsg <- function(msg) {
    cat(msg, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  # Log input parameters
  logMsg("Input parameters:\n")
  logMsg(paste0("  Reference raster: ", referenceRasterFile, "\n"))
  logMsg("  Explanatory variable rasters:\n")
  for (i in seq_along(varRasterFiles)) { logMsg(paste0("  [", i, "] ", varRasterFiles[[i]], "\n")) }
  logMsg(paste0("  Initiation points: ", initiationPointsFile, "\n"))
  logMsg(paste0("  Non-initiation points ratio: ", noninitiationRatio, "\n"))
  logMsg(paste0("  Buffer radius: ", bufferRadius, "\n"))
  logMsg(paste0("  Buffer extraction method: ", bufferExtractionMethod, "\n"))
  logMsg(paste0("  Initiation limit percent: ", initiationLimitPercent, "%\n"))
  logMsg(paste0("  k: ", k, "\n"))
  logMsg(paste0("  Generate probability rasters: ", generateProbabilityRasters, "\n"))
  logMsg(paste0("  Output directory: ", outputDir, "\n\n"))
  
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
  initiationRange <- createInitiationRange(varsRaster, initiationBuffers, initiationLimitPercent)
  
  # The region where non-initiation buffers can be sampled from: regions that 
  # meet initiation conditions but recorded no landslides
  noninitiationRaster <- createNoninitiationRaster(
    varsRaster,
    initiationRange,
    initiationBuffers
  )
  
  # Generate non-initiation buffers
  noninitiationBuffersCount <- ceiling(length(initiationPoints) * noninitiationRatio)
  noninitiationBuffers <- generateNoninitiationBuffers(
    noninitiationBuffersCount,
    noninitiationRaster,
    bufferRadius
  )
  
  # Plot variable distributions ------------------------------------------------
  
  # Get requested values from all initiation and non-initiation buffers
  allBuffersData <- extractBufferValues(
    varsRaster,
    initiationBuffers,
    noninitiationBuffers,
    bufferExtractionMethod
  )
  
  # For each explanatory variable, draw a box plot of all the buffer values
  for (varName in names(varsRaster)) {
    varInitiationValues <- allBuffersData[allBuffersData$class == "initiation", varName]
    varNoninitiationValues <- allBuffersData[allBuffersData$class == "non-initiation", varName]
    
    # Save plot as an image file
    dev.new()
    boxplot(
      varInitiationValues, varNoninitiationValues,
      main = paste0(varName, " distribution"),
      names = c("Initiation", "Non-initiation"),
      col = c("green", "red")
    )
    dev.copy(jpeg, paste0(outputDir, "/", varName, "_dist.jpeg"))
    dev.off()
  }
  
  # Create initiation buffer sets ----------------------------------------------
  
  # Calculate how many testing initiation buffers there should be per iteration
  testingInitiationBuffersPerIteration <- floor(length(initiationBuffers) / k)
  
  # Create k sets of initiation buffers to use for training and testing
  testingSets <- list()
  testingFreeInitiationIndices <- seq_along(initiationBuffers)
  for (i in 1:k) {
    testingInitiationIndices <- sample(testingFreeInitiationIndices, size = testingInitiationBuffersPerIteration)
    testingFreeInitiationIndices <- testingFreeInitiationIndices[!(testingFreeInitiationIndices %in% testingInitiationIndices)]
    testingSets[[i]] <- testingInitiationIndices
  }
  
  # Create non-initiation buffer sets ------------------------------------------
  
  # Create a static set of training and testing non-initiation buffers to use every iteration
  testingNoninitiationCount <- floor(testingInitiationBuffersPerIteration * noninitiationRatio)
  testingNoninitiationIndices <- sample(seq_along(noninitiationBuffers), size = testingNoninitiationCount)
  trainingNoninitiationIndices <- setdiff(seq_along(noninitiationBuffers), testingNoninitiationIndices)
  
  trainingNoninitiationBuffers <- noninitiationBuffers[trainingNoninitiationIndices]
  testingNoninitiationBuffers <- noninitiationBuffers[testingNoninitiationIndices]
  
  # Perform cross-validation ---------------------------------------------------
  
  # Place to store iteration model auc values
  iterationsAucValues <- c()
  
  # Place to store iteration model error rates
  iterationsErrorRates <- data.frame(
    rep(NA, k),
    rep(NA, k),
    rep(NA, k)
  )
  names(iterationsErrorRates) <- c("OOB", "initiation", "non-initiation")
  
  for (i in seq_along(testingSets)) {
    ## Create training data ----------------------------------------------------
    
    testingInitiationIndices <- testingSets[[i]]
    
    # Get training initiation buffers
    trainingInitiationIndices <- setdiff(seq_along(initiationBuffers), testingInitiationIndices)
    trainingInitiationBuffers <- initiationBuffers[trainingInitiationIndices]
    
    # Create training dataset
    trainingData <- extractBufferValues(
      varsRaster,
      trainingInitiationBuffers,
      trainingNoninitiationBuffers,
      bufferExtractionMethod
    )
    
    ## Create model ------------------------------------------------------------
    
    # Train a new random forest model
    rfModel <- randomForest::randomForest(
      formula = class ~ .,
      data = trainingData
    )
    
    logMsg(paste0("Model ", formatC(i, width = 2, format = "d", flag = "0"), 
                  " -------------------------------------------\n\n"))
    
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
    
    ## Run model on testing data -----------------------------------------------
    
    # Get testing initiation buffers
    testingInitiationBuffers <- initiationBuffers[testingInitiationIndices]
    
    # Creating testing data
    testingData <- extractBufferValues(
      varsRaster,
      testingInitiationBuffers,
      testingNoninitiationBuffers,
      bufferExtractionMethod
    )
    
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
    
    # Record iteration statistics
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
  colnames(aucMatrix) <- c("min", "max", "range", "stdev")
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
      noninitiationRatio = 1.5,
      bufferRadius = 20,
      bufferExtractionMethod = "all cells",
      initiationLimitPercent = 10,
      k = 5,
      generateProbabilityRasters = FALSE,
      outputDir = "C:/Work/netmapdata/Scottsburg"
    ),
    out_params = list()
  )
  
  # Run in Scottsburg (DESKTOP2)
  tool_exec(
    in_params = list(
      referenceRasterFile = "E:/NetmapData/Scottsburg/elev_scottsburg.flt",
      varRasterFiles = list(
        "E:/NetmapData/Scottsburg/grad_30.tif",
        "E:/NetmapData/Scottsburg/plan_30.tif"
      ),
      initiationPointsFile = "E:/NetmapData/Scottsburg/Scottsburg_Upslope.shp",
      noninitiationRatio = 1.5,
      bufferRadius = 20,
      bufferExtractionMethod = "all cells",
      initiationLimitPercent = 10,
      k = 5,
      generateProbabilityRasters = FALSE,
      outputDir = "E:/NetmapData/Scottsburg"
    ),
    out_params = list()
  )
  
 }
 