#' @title Perform cross-validation
#'
#' @description Performs repeated k-fold spatial cross-validation using 
#' different sets of generated non-initiation points. Outputs an average 
#' landslide susceptibility raster and a log file that records and summarizes 
#' model performance.
#' 
#' The process:
#' 1. Use initiation buffers to define an analysis region
#' 2. Generate a set of non-initiation buffers in the analysis region
#' 3. Extract landslide dataset from all buffers
#' 4. Perform repeated k-fold spatial cross-validation on dataset
#' 5. Repeat from step 2 a given number of times 
#'
#' @param refRasterFile          A raster file to use as a grid reference.
#' @param varRasterFiles         A list of raster files to use as explanatory 
#'                               variables.
#' @param initPointsFile         A shapefile of initiation points.
#' @param noninitRatio           The ratio of non-initiation sites to initiation
#'                               sites.
#' @param bufferRadius           The radius of site buffers.
#' @param bufferExtractionMethod The method used to select values from site 
#'                               buffers for training/testing. Either: 
#'                               "all cells", "center cell", 
#'                               "max gradient cell", or "max plan cell".
#' @param initRangeExpansion     The proportion (in %) to expand the initiation 
#'                               range of each variable by.
#' @param noninitSetsCount       The number of different non-initiation sets to
#'                               generate.
#' @param repetitionsCount       The number of times to repeat k-fold cross-
#'                               validation.
#' @param foldsCount             The number of folds to use for k-fold cross-
#'                               validation.
#' @param generateAvgProbRaster  Should an average probability raster be 
#'                               generated?
#' @runName                      Name of this run. Used to name generated output
#'                               files.
#' @param outputDir              The directory to write output files to.
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
#'   noninitSetsCount       = 1, 
#'   repetitionsCount       = 1,
#'   foldsCount             = 4,
#'   generateAvgProbRaster  = FALSE,
#'   runName                = "scott_cv",
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
  repetitionsCount,
  foldsCount,
  generateAvgProbRaster,
  runName,
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
  source("./helper/summarizeVector.R")
  source("./helper/summarizeDataFrame.R")
  
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
    stop("Non-initiation points ratio cannot be greater than 0.")
  
  # Validate buffer radius
  if (bufferRadius < 0)
    stop("Buffer radius must be greatr than or equal to 0.")
  
  # Validate buffer extraction method
  if (!(bufferExtractionMethod %in% c("all cells", "center cell", "max gradient 
      cell", "max plan cell")))
    stop("Buffer extraction method must be one of 'all cells', 'center cell', 
       'max gradient cell', or 'max plan cell'.")
  
  # Validate number of non-initiation sets
  if (noninitSetsCount < 1)
    stop("Non-initiation sets must be greater than or equal to 1.")
  
  # Validate number of repetitions
  if (repetitionsCount < 1)
    stop("Repetitions must be greater than or equal to 1.")
  
  # Validate number of folds
  if (foldsCount < 2)
    stop("Folds must be greater than 1.")
  
  if (!is.logical(generateAvgProbRaster))
    generateAvgProbRaster <- FALSE
  
  # Validate run name
  if (nchar(runName) < 1)
    stop("Must provide a run name")
  
  # Validate output directory
  if (!file.exists(outputDir))
    stop("Output directory not found: '", outputDir, "'.")
  
  # Set up logging -------------------------------------------------------------
  
  logFile <- paste0(outputDir, "/", runName, "_log.txt")
  file.create(logFile)
  
  logMsg <- function(msg) cat(msg, file = logFile, append = TRUE)
  logObj <- function(obj) capture.output(obj, file = logFile, append = TRUE)
  
  # Log input parameters
  logMsg("INPUT PARAMETERS:\n")
  logMsg(paste0("  Reference raster: ", refRasterFile, "\n"))
  logMsg("  Explanatory variable rasters:\n")
  for (i in seq_along(varRasterFiles)) 
    logMsg(paste0("    [", i, "] ", varRasterFiles[[i]], "\n"))
  logMsg(paste0("  Initiation points: ", initPointsFile, "\n"))
  logMsg(paste0("  Non-initiation points ratio: ", noninitRatio, "\n"))
  logMsg(paste0("  Buffer radius: ", bufferRadius, "\n"))
  logMsg(paste0("  Buffer extraction method: ", bufferExtractionMethod, "\n"))
  logMsg(paste0("  Initiation range expansion: ", initRangeExpansion, "%\n"))
  logMsg(paste0("  # non-initiation points sets: ", noninitSetsCount, "\n"))
  logMsg(paste0("  # k-fold CV repetitions: ", repetitionsCount, "\n"))
  logMsg(paste0("  # folds (k): ", foldsCount, "\n"))
  logMsg(paste0("  Generate probability raster? ", generateAvgProbRaster, "\n"))
  logMsg(paste0("  Run name: ", runName, "\n"))
  logMsg(paste0("  Output directory: ", outputDir, "\n"))
  logMsg("\n")
  
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
  
  if (generateAvgProbRaster) {
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
    
    # Remove coordinates columns from landslide dataset
    coordsCols <- names(landslideData) %in% c("x", "y")  
    landslideData <- landslideData[,!coordsCols]
    
    # Perform repeated k-fold spatial cross-validation -------------------------
    
    for (repetition in seq_len(resampling$iters)) {
      
      ## Train model -----------------------------------------------------------
      
      # Get training data
      trainingIndices <- resampling$train_set(repetition)
      trainingData <- landslideData[trainingIndices,]
      
      # Train random forest model
      rfModel <- randomForest::randomForest(
        formula = class ~ .,
        data = trainingData
      )
      
      logMsg(paste0("Model ", formatC(iterationNumber, width = 3, format = "d", 
        flag = "0"), " ------------------------------------------------\n\n"))
      
      # Log model error rates
      logMsg("TRAINING ERROR RATES:\n")
      trainingEr <- data.frame(rfModel$err.rate[rfModel$ntree,])
      colnames(trainingEr) <- "error rate"
      logObj(trainingEr)
      logMsg("\n")
      
      ## Test model ------------------------------------------------------------
      
      # Get testing data
      testingIndices <- resampling$test_set(repetition)
      testingData <- landslideData[testingIndices,]
      
      # Have model classify testing data
      prediction <- predict(
        rfModel,
        type = "response",
        newdata = testingData
      )
      
      # Have model predict initiation probability of test data
      initProb <- predict(
        rfModel,
        type = "prob",
        newdata = testingData
      )[,"initiation"]
      
      ## Record model test results ---------------------------------------------
      
      # Calculate ROC stats
      rocStats <- TerrainWorksUtils::calcRocStats(
        classes = testingData$class,
        probs = initProb,
        "initiation",
        "non-initiation"
      )
      
      aucValue <- rocStats$auc@y.values[[1]]
      iterationsAucValues <- c(iterationsAucValues, aucValue)
      
      errorRates <- as.data.frame(rfModel$err.rate)[rfModel$ntree,]
      iterationsErrorRates <- rbind(iterationsErrorRates, errorRates)
      
      # Log test confusion matrix
      logMsg("TESTING CONFUSION MATRIX:\n")
      testCf <- table(prediction, testingData$class)
      testCf <- data.frame(testCf[1,], testCf[2,])
      colnames(testCf) <- c("initiation", "non-initiation")
      rownames(testCf) <- c("true initiation", "true non-initiation")
      logObj(testCf)
      logMsg("\n")
      
      
      logMsg(paste0("TESTING AUC: ", round(aucValue, digits = 7), "\n"))
      logMsg("\n")
      
      ## Generate probability raster -------------------------------------------
      
      if (generateAvgProbRaster) {
        # Have the model predict an initiation probability raster
        probRaster <- terra::predict(
          analysisRegionVarsRaster,
          rfModel,
          na.rm = TRUE,
          type = "prob"
        )[["initiation"]]
        
        # Save the probability raster to disk
        probRasterFile <- tempfile("prob", outputDir, ".tif")
        terra::writeRaster(probRaster, probRasterFile)
        
        # Record the name of the probability raster file
        probRasterFiles[iterationNumber] <- probRasterFile
      }
      
      iterationNumber <- iterationNumber + 1
      
    }
    
  }
  
  # Summarize test results -----------------------------------------------------
  
  logMsg("SUMMARY --------------------------------------------------\n\n")
  
  # Summarize error rates
  errorRatesSummary <- summarizeDataFrame(iterationsErrorRates)
  
  # Summarize AUC values
  aucValuesSummary <- summarizeVector(iterationsAucValues)
  rownames(aucValuesSummary) <- "AUC"
  
  # Log summary of error rates
  logMsg("TESTING ERROR RATES:\n")
  logObj(errorRatesSummary)
  logMsg("\n")
  
  # Log summary of AUC values
  logMsg("TESTING AUC:\n")
  logObj(aucValuesSummary)
  logMsg("\n")
  
  # Generate average probability raster ----------------------------------------
  
  if (generateAvgProbRaster) {
    # Create average probability raster
    avgProbRasterFile <- paste0(outputDir, "/", runName, "_prob.tif")
    avgProbRaster <- createAverageRaster(probRasterFiles, avgProbRasterFile)
    
    # Save average probability raster
    terra::writeRaster(avgProbRaster, avgProbRasterFile, overwrite = TRUE)
    
    # Delete temporary probability raster files
    sapply(probRasterFiles, function(file) unlink(file))
  }
  
}

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
    repetitionsCount       = in_params[[9]],
    foldsCount             = in_params[[10]],
    generateAvgProbRaster  = in_params[[11]], 
    runName                = in_params[[12]],
    outputDir              = in_params[[13]]
  )
  
  return(out_params)
  
}
