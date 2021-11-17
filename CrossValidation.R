#' @title Perform cross-validation
#'
#' @description Performs repeated k-fold spatial cross-validation using 
#' different sets of generated non-initiation points. Outputs a log file that 
#' records and summarizes model performance.
#' 
#' The process:
#' 1. Use initiation buffers to define an analysis region
#' 2. Generate a set of non-initiation buffers in the analysis region
#' 3. Extract landslide data values from all buffers
#' 4. Perform repeated k-fold spatial cross-validation on dataset
#' 5. Repeat from step 2 a given number of times 
#'
#' @param refRasterFile          File of a raster to use as a grid reference 
#' @param varRasterFiles         List of raster files to use as explanatory 
#'                               variables
#' @param initPointsFile         File of initiation points
#' @param noninitRatio           The ratio of non-initiation sites to initiation
#'                               sites
#' @param bufferRadius           The radius of site buffers
#' @param bufferExtractionMethod The method to select values from site buffers. 
#'                               Either: "all cells", "center cell", "max 
#'                               gradient cell", or "max plan cell"
#' @param initRangeExpansion     The proportion (in %) to expand each variable 
#'                               initiation range by
#' @param noninitSetsCount       The number of different non-initiation sets to
#'                               generate
#' @param foldsCount             The number of folds to use for k-fold cross-
#'                               validation
#' @param repetitionsCount       How many times to repeat k-fold cross-
#'                               validation for each non-initiation set
#' @param generateProbRaster     Should an average probability raster be 
#'                               generated?
#' @param outputDir              The directory to write output files to
#' 
#' @example
#' \donttest {
#' performCrossValidation(
#'   refRasterFile          = "E:/NetmapData/Scottsburg/elev_scottsburg.flt",
#'   varRasterFiles         = list(
#'     "E:/NetmapData/Scottsburg/grad_30.tif",
#'     "E:/NetmapData/Scottsburg/plan_30.tif"
#'   ),
#'   initPointsFile         = "E:/NetmapData/Scottsburg/Scottsburg_Upslope.shp",
#'   noninitRatio           = 1,
#'   bufferRadius           = 20,
#'   bufferExtractionMethod = "center cell",
#'   initRangeExpansion     = 0,
#'   noninitSetsCount       = 10, 
#'   foldsCount             = 4,
#'   repetitionsCount       = 5,
#'   generateProbRaster     = true,
#'   outputDir              = "E:/NetmapData/Scottsburg"
#' )
#' }

performCrossValidation <- function(
  refRasterFile,
  varRasterFiles,
  initPointsFile,
  noninitRatio,
  bufferRadius,
  bufferExtractionMethod,
  initRangeExpansion,
  noninitSetsCount,
  foldsCount,
  repetitionsCount,
  generateProbRaster,
  outputDir
) {
  
  # Install dependencies -------------------------------------------------------
  
  dependencies <- c(
    "mlr3",
    "mlr3spatiotempcv",
    "terra"
  )
  
  install.packages(setdiff(dependencies, rownames(installed.packages())))
  
  # Load helper functions ------------------------------------------------------
  
  source("./helper/createInitiationRange.R")
  source("./helper/createAnalysisRegionMask.R")
  source("./helper/generateNoninitiationBuffers.R")
  source("./helper/extractBufferValues.R")
  
  # Validate parameters --------------------------------------------------------
  
  # Validate reference raster file
  if (!file.exists(refRasterFile))
    stop(paste0("Reference raster file not found: '", refRasterFile, "'."))
  
  # Validate variable raster files
  lapply(varRasterFiles, function(file) {
    if (!file.exists(file))
      stop(paste0("Variable raster file not found: '", file, "'."))
  })
  
  # Validate initiation points file
  if (!file.exists(initPointsFile))
    stop(paste0("Initiation points file not found: '", initPointsFile, "'."))
  
  # Validate non-initiation ratio
  if (noninitRatio <= 0)
    stop("Non-initiation points ratio cannot be <=0.")
  
  # Validate buffer radius
  if (bufferRadius < 0)
    stop("Buffer radius cannot be negative.")
  
  # Validate buffer extraction method
  if (!(bufferExtractionMethod %in% c("all cells", "center cell", "max gradient 
      cell", "max plan cell")))
    stop("Buffer extraction method must be one of 'all cells', 'center cell', 
       'max gradient cell', or 'max plan cell'.")
  
  # Validate number of non-initiation sets
  if (noninitSetsCount < 1)
    stop("Non-initiation sets cannot be fewer than 1.")
  
  # Validate number of folds
  if (foldsCount < 1)
    stop("Folds cannot be fewer than 1.")
  
  # Validate number of repetitions
  if (repetitionsCount < 1)
    stop("Repetitions cannot be fewer than 1.")
  
  # Validate output directory
  if (!file.exists(outputDir))
    stop("Output directory not found: '", outputDir, "'.")
  
  # Set up logging -------------------------------------------------------------
  
  logFilename <- "cv_log.txt"
  file.create(paste0(outputDir, "/", logFilename))
  
  logMsg <- function(msg) {
    cat(msg, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  logObj <- function(obj) {
    capture.output(obj, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  # Load rasters ---------------------------------------------------------------
  
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
  initRange <- createInitiationRange(varsRaster, initBuffers, initRangeExpansion)
  
  # Log initiation range matrix
  logMsg("INITIATION RANGES:\n")
  logObj(initRange)
  logMsg("\n")
  
  # Identify cells in the study region that have variable values within their 
  # initiation ranges
  analysisRegionMask <- createAnalysisRegionMask(varsRaster, initRange)
  
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
  
  # Perform cross-validation using different sets of non-initiation buffers ----
  
  # Total number of iterations
  iterationsCount <- noninitSetsCount * foldsCount * repetitionsCount
  
  # Place to store model auc value for each iteration
  iterationsAucValues <- c()
  
  # Place to store model error rates for each iteration
  iterationsErrorRates <- data.frame(
    rep(NA, iterationsCount),
    rep(NA, iterationsCount),
    rep(NA, iterationsCount)
  )
  names(iterationsErrorRates) <- c("OOB", "initiation", "non-initiation")
  
  if (generateProbRaster) {
    # Create a version of the variables raster that only keeps cells within the 
    # analysis region
    analysisRegionVarsRaster <- terra::mask(varsRaster, analysisRegionMask)
    
    # Vector of temporary probability raster files
    probRasterFiles <- c()
  }
  
  iterationNumber <- 1
  
  for (noninitSetNumber in seq_len(noninitSetsCount)) {
    
    # Generate non-initiation buffers ------------------------------------------
    
    noninitBuffers <- generateNoninitiationBuffers(
      noninitBuffersCount,
      noninitRegion,
      bufferRadius
    )
    
    # Create full landslide dataset --------------------------------------------
    
    landslideData <- extractBufferValues(
      varsRaster,
      initBuffers,
      noninitBuffers,
      bufferExtractionMethod
    )
    
    # Define cross-validation method -------------------------------------------
    
    # Set up a machine learning task for landslide data
    task <- mlr3spatiotempcv::TaskClassifST$new(
      "landslides",
      backend = mlr3::as_data_backend(landslideData),
      target = "class",
      positive = "initiation",
      extra_args = list(
        coordinate_names = c("x", "y"),
        crs = terra::crs(refRaster, proj = TRUE)
      )
    )
    
    # Create a version of the variables raster that only keeps cells within the 
    # analysis region
    analysisVarsRaster <- terra::mask(varsRaster, analysisRegionMask)
    
    # Set up a resampling method for repeated k-fold spatial cross-validation
    resampling <- mlr3::rsmp(
      "repeated_spcv_coords",
      folds = foldsCount,
      repeats = repetitionsCount
    )
    
    # Perform the resampling method on the task
    resampling <- resampling$instantiate(task)
    
    # Perform repeated k-fold spatial cross-validation -------------------------
    
    for (repetition in seq_len(resampling$iters)) {
      
      ## Train model -----------------------------------------------------------
      
      # Get training data
      trainingIndices <- resampling$train_set(repetition)
      trainingData <- landslideData[trainingIndices,]
      
      # Remove columns for coordinates
      coordsCols <- names(trainingData) %in% c("x", "y")  
      trainingData <- trainingData[,!coordsCols]
      
      # Train random forest model
      rfModel <- randomForest::randomForest(
        formula = class ~ .,
        data = trainingData
      )
      
      logMsg(paste0("Model ", formatC(iterationNumber, width = 3, format = "d", 
        flag = "0"), " ------------------------------------------------\n\n"))
      
      # Log model error rates
      logMsg("MODEL ERROR RATES:\n")
      errorRateDf <- data.frame(rfModel$err.rate[rfModel$ntree,])
      colnames(errorRateDf) <- "error rate"
      logObj(errorRateDf)
      logMsg("\n")
      
      ## Evaluate model --------------------------------------------------------
      
      # Get testing data
      testingIndices <- resampling$test_set(repetition)
      testingData <- landslideData[testingIndices,]
      
      # Remove columns for coordinates
      coordsCols <- names(testingData) %in% c("x", "y")
      testingData <- testingData[,!coordsCols]
      
      # # Predict the class type of each test dataset entry
      prediction <- predict(
        rfModel,
        type = "response",
        newdata = testingData
      )

      # Log test dataset confusion matrix
      logMsg("TEST RESULTS:\n")
      testConfusionMatrix <- table(prediction, testingData$class)
      testConfusionMatrix <- data.frame(testConfusionMatrix[1,], testConfusionMatrix[2,])
      colnames(testConfusionMatrix) <- c("initiation", "non-initiation")
      rownames(testConfusionMatrix) <- c("true initiation", "true non-initiation")
      logObj(testConfusionMatrix)
      logMsg("\n")
      
      # Calculate initiation probability for each test entry
      initProb <- predict(
        rfModel,
        type = "prob",
        newdata = testingData
      )[,"initiation"]
      
      # Calculate ROC stats
      rocStats <- TerrainWorksUtils::calcRocStats(
        classes = testingData$class,
        probs = initProb,
        "initiation",
        "non-initiation"
      )
      
      ## Record iteration statistics -------------------------------------------
      
      aucValue <- rocStats$auc@y.values[[1]]
      iterationsAucValues <- c(iterationsAucValues, aucValue)
      
      errorRates <- as.data.frame(rfModel$err.rate)[rfModel$ntree,]
      iterationsErrorRates <- rbind(iterationsErrorRates, errorRates)
      
      logMsg(paste0("AUC: ", round(aucValue, digits = 7), "\n"))
      logMsg("\n")
      
      ## Generate probability raster -------------------------------------------
      
      if (generateProbRaster) {
        # Generate a probability raster using this model
        probRaster <- terra::predict(
          analysisRegionVarsRaster,
          rfModel,
          na.rm = TRUE,
          type = "prob"
        )[["initiation"]]
        
        # Save the probability raster to disk
        probRasterFile <- tempfile("prob", outputDir, ".tif")
        terra::writeRaster(probRaster, probRasterFile)
        
        # Record the name of the probability raster
        probRasterFiles[iterationNumber] <- probRasterFile 
      }
      
      ## Increment the iteration -----------------------------------------------
      
      iterationNumber <- iterationNumber + 1
      
    }
    
  }
  
  # Summarize model performance ------------------------------------------------
  
  logMsg("SUMMARY --------------------------------------------------\n\n")
  
  # Summarize AUC value and error rates across all iterations
  aucValuesSummary <- summarizeAucValues(iterationsAucValues)
  errorRatesSummary <- summarizeErrorRates(iterationsErrorRates)
  
  # Log summary of AUC values
  logMsg("AUC:\n")
  logObj(aucValuesSummary)
  logMsg("\n")
  
  # Log summary of error rates
  logMsg("ERROR RATES:\n")
  logObj(errorRatesSummary)
  
  # Generate average probability raster ----------------------------------------
  
  if (generateProbRaster) {
    # Create average probability raster
    avgProbRasterFile <- paste0(outputDir, "/avg_prob.tif")
    avgProbRaster <- createAverageRaster(probRasterFiles, avgProbRasterFile)
    
    # Save average probability raster
    terra::writeRaster(avgProbRaster, avgProbRasterFile, overwrite = TRUE)
    
    # Delete temporary probability raster files
    sapply(probRasterFiles, function(file) unlink(file))
  }
  
}

summarizeAucValues <- function(aucValues) {
  
  # Calculate AUC value statistics
  aucMin <- min(aucValues, na.rm = TRUE)
  aucMax <- max(aucValues, na.rm = TRUE)
  aucRange <- diff(c(aucMin, aucMax))
  aucMean <- mean(aucValues, na.rm = TRUE)
  aucSd <- sd(aucValues, na.rm = TRUE)
  
  # Create summary table
  aucSummary <- matrix(c(aucMin, aucMax, aucRange, aucMean, aucSd), ncol = 1)
  aucSummary <- t(aucSummary)
  colnames(aucSummary) <- c("min", "max", "range", "mean", "sd")
  rownames(aucSummary) <- "AUC"
  
  return(aucSummary)
  
}

summarizeErrorRates <- function(errorRates) {
  
  # Calculate error rates statistics
  erMin <- as.data.frame(lapply(errorRates, min, na.rm = TRUE))
  erMax <- as.data.frame(lapply(errorRates, max, na.rm = TRUE))
  erRange <- as.data.frame(lapply(errorRates, function(x) {
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  }))
  erMean <- as.data.frame(lapply(errorRates, mean, na.rm = TRUE))
  erSd <- as.data.frame(lapply(errorRates, sd, na.rm = TRUE))
  
  # Create summary table
  erDf <- rbind(erMin, erMax, erRange, erMean, erSd)
  erDf <- t(erDf)
  erSummary <- as.matrix(erDf)
  colnames(erSummary) <- c("min", "max", "range", "mean", "sd")
  
  return(erSummary)
  
}

# Averages a set of rasters on disk and saves the result. Necessary
# when low on memory or default temporary disk space.
createAverageRaster <- function(rasterFiles, filename) {
  
  if (length(rasterFiles) == 0)
    return(NULL)
  
  sumRaster <- terra::rast(rasterFiles[1])
  if (length(rasterFiles) > 1) {
    for (i in 2:length(rasterFiles)) {
      sumRaster <- sumRaster + terra::rast(rasterFiles[i])
    }
  }
  avgRaster <- sumRaster / length(rasterFiles)
    
  return(avgRaster)
  
}

# Entrypoint for ArcGIS script tool
tool_exec <- function(in_params, out_params) {
  
  performCrossValidation(
    refRasterFile          = in_params[[1]],
    varRasterFiles         = in_params[[2]],
    initPointsFile         = in_params[[3]],
    noninitRatio           = in_params[[4]],
    bufferRadius           = in_params[[5]],
    bufferExtractionMethod = in_params[[6]],
    initRangeExpansion     = in_params[[7]],
    noninitSetsCount       = in_params[[8]],
    foldsCount             = in_params[[9]],
    repetitionsCount       = in_params[[10]],
    generateProbRaster     = in_params[[11]],
    outputDir              = in_params[[12]]
  )
  
  return(out_params)
  
}
