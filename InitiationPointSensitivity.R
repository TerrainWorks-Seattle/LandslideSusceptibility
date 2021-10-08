# This script provides a way to assess how much the choice of initiation points
# affects the accuracy of a landslide susceptibility prediction model.
# k iterations of the model are created, each trained with a unique set of
# initiation points and a static set of non-initiation points. A buffer of given
# size is drawn around each point and all cell values within are used as model 
# inputs.

tool_exec <- function(in_params, out_params) {
  
  # Load helper functions ------------------------------------------------------
  
  source("./createInitiationRange.R")
  source("./createNoninitiationRaster.R")
  source("./generateNoninitiationBuffers.R")
  source("./extractBufferValues.R")
  
  # Set parameters -------------------------------------------------------------
  
  # Input
  referenceRasterFile        <- in_params[[1]]
  varRasterFiles             <- in_params[[2]]
  initiationPointsFile       <- in_params[[3]]
  noninitiationRatio         <- in_params[[4]]
  bufferRadius               <- in_params[[5]]
  bufferExtractionMethod     <- in_params[[6]]
  initiationRangeExpansion   <- in_params[[7]]
  iterations                 <- in_params[[8]]
  testingProportionPercent   <- in_params[[9]]
  generateProbabilityRasters <- in_params[[10]]
  outputDir                  <- in_params[[11]]
  
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
  logMsg(paste0("  Testing proportion: ", testingProportionPercent, "%\n"))
  logMsg(paste0("  Generate probability rasters: ", generateProbabilityRasters, "\n"))
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
  
  # Create non-initiation buffer sets ------------------------------------------
  
  # Create a static set of training and testing non-initiation buffers to use every iteration
  testingNoninitiationCount <- floor(length(noninitiationBuffers) * (testingProportionPercent / 100))
  testingNoninitiationIndices <- sample(seq_along(noninitiationBuffers), size = testingNoninitiationCount)
  trainingNoninitiationIndices <- setdiff(seq_along(noninitiationBuffers), testingNoninitiationIndices)
  
  trainingNoninitiationBuffers <- noninitiationBuffers[trainingNoninitiationIndices]
  testingNoninitiationBuffers <- noninitiationBuffers[testingNoninitiationIndices]
  
  # Perform cross-validation ---------------------------------------------------
  
  # Place to store iteration model AUC values
  iterationsAucValues <- c()
  
  # Place to store iteration model error rates
  iterationsErrorRates <- data.frame(
    rep(NA, iterations),
    rep(NA, iterations),
    rep(NA, iterations)
  )
  names(iterationsErrorRates) <- c("OOB", "initiation", "non-initiation")
  
  # Calculate how many initiation buffers should be used for testing per iteration
  testingInitiationBuffersCount <- floor(length(initiationBuffers) * (testingProportionPercent / 100))
  
  for (i in seq_len(iterations)) {
    ## Create initiation buffer sets -------------------------------------------
    
    # Create testing initiation set
    testingInitiationIndices <- sample(seq_along(initiationBuffers), size = testingInitiationBuffersCount)
    testingInitiationBuffers <- initiationBuffers[testingInitiationIndices]
    
    # Create training initiation set
    trainingInitiationIndices <- setdiff(seq_along(initiationBuffers), testingInitiationIndices)
    trainingInitiationBuffers <- initiationBuffers[trainingInitiationIndices]
    
    ## Create model ------------------------------------------------------------
    
    # Create training dataset
    trainingData <- extractBufferValues(
      varsRaster,
      trainingInitiationBuffers,
      trainingNoninitiationBuffers,
      bufferExtractionMethod
    )
    
    # Train a new random forest model
    rfModel <- randomForest::randomForest(
      formula = class ~ .,
      data = trainingData
    )
    
    logMsg(paste0("Model ", formatC(i, width = 2, format = "d", flag = "0"), 
                  " ------------------------------------------------\n\n"))
    
    # Log model error rates
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
    
    ## Run model on testing data -----------------------------------------------
    
    # Creating testing data
    testingData <- extractBufferValues(
      varsRaster,
      testingInitiationBuffers,
      testingNoninitiationBuffers,
      bufferExtractionMethod
    )
    
    # Predict the class type of each test dataset entry
    prediction <- predict(
      rfModel,
      type = "response",
      newdata = testingData
    )
    
    # Log test dataset confusion matrix
    logMsg("TESTING CONFUSION MATRIX:\n")
    testConfusionMatrix <- table(prediction, testingData$class)
    testConfusionMatrix <- data.frame(testConfusionMatrix[1,], testConfusionMatrix[2,])
    colnames(testConfusionMatrix) <- c("initiation", "non-initiation")
    rownames(testConfusionMatrix) <- c("true initiation", "true non-initiation")
    logObj(testConfusionMatrix)
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
    logMsg(paste0("TESTING AUC: ", round(auc, digits = 7), "\n"))
    logMsg("\n")
    
    # Record iteration statistics
    iterationsAucValues <- c(iterationsAucValues, auc)
    iterationsErrorRates[i,] <- as.data.frame(rfModel$err.rate)[rfModel$ntree,]
    
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
      bufferExtractionMethod = "center cell",
      initiationRangeExtension = 10,
      k = 5,
      testingProportionPercent = 10,
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
      bufferExtractionMethod = "max gradient cell",
      initiationRangeExpansion = 10,
      k = 20,
      testingProportionPercent = 10,
      generateProbabilityRasters = FALSE,
      outputDir = "E:/NetmapData/Scottsburg"
    ),
    out_params = list()
  )
  
 }
 