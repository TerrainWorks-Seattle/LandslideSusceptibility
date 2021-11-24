#' @title Assess the sensitivity of landslide susceptibility models to the 
#' choice of initiation sites
#'
#' @description Creates a series of landslide susceptibility models trained on 
#' random subsets of initiation sites. Outputs a distribution plot of site 
#' values for each variable and a log file that records model accuracy. 
#'
#' @param refRasterFile          File of a raster to use as a grid reference 
#' @param varsRasterFile         List of raster files to use as explanatory 
#'                               variables
#' @param initPointsFile         File of initiation points
#' @param noninitRatio           The ratio of non-initiation sites to initiation
#'                               sites
#' @param bufferRadius           The radius of site buffers
#' @param bufferExtractionMethod The method to select values from site buffers. 
#'                               Either: "all cells", "center cell", "max 
#'                               gradient cell", or "max plan cell"
#' @param initRangeExpansion     The relative proportion (in %) to expand each 
#'                               variable initiation range by
#' @param iterations             How many models to create
#' @param testingProportion      The proportion (in %) of initiation and
#'                               non-initiation sites to withhold for testing
#' @param outputDir              The directory to write output files to
#' 
#' @example
#' \donttest{
#' assessInitiationPointSusceptibility(
#'   refRasterFile = "E:/NetmapData/Scottsburg/elev_scottsburg.flt",
#'   varRasterFiles = list(
#'     "E:/NetmapData/Scottsburg/grad_30.tif",
#'     "E:/NetmapData/Scottsburg/plan_30.tif",
#'     "E:/NetmapData/Scottsburg/pca_scott.flt"
#'   ),
#'   initPointsFile = "E:/NetmapData/Scottsburg/Scottsburg_Upslope.shp",
#'   noninitRatio = 1.5,
#'   bufferRadius = 20,
#'   bufferExtractionMethod = "max gradient cell",
#'   initRangeExpansion = 10,
#'   iterations = 20,
#'   testingProportion = 10,
#'   outputDir = "E:/NetmapData/Scottsburg"
#' )
#' }

assessInitiationPointSusceptibility <- function(
  refRasterFile,
  varRasterFiles,
  initPointsFile,
  noninitRatio,
  bufferRadius,
  bufferExtractionMethod,
  initRangeExpansion,
  iterationsCount,
  testingProportion,
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
    stop("Non-initiation points ratio cannot be <=0.")
  
  # Validate buffer radius
  if (bufferRadius < 0)
    stop("Buffer radius cannot be negative.")
  
  # Validate buffer extraction method
  if (!(bufferExtractionMethod %in% c("all cells", "center cell", 
                                      "max gradient cell", "max plan cell")))
    stop("Buffer extraction method must be one of 'all cells', 'center cell', 
         'max gradient cell', or 'max plan cell'.")
  
  # Validate number of iterations
  if (iterationsCount < 1)
    stop("Iterations cannot be fewer than 1.")
  
  # Validate testing proportion
  if (testingProportion <= 0 || testingProportion >= 100)
    stop("Testing proportion must be between 0% and 100%.")
  
  # Validate output directory
  if (!file.exists(outputDir))
    stop("Output directory not found: '", outputDir, "'.")
  
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
  logMsg(paste0("  Reference raster: ", refRasterFile, "\n"))
  logMsg("  Explanatory variable rasters:\n")
  for (i in seq_along(varRasterFiles)) 
    logMsg(paste0("    [", i, "] ", varRasterFiles[[i]], "\n"))
  logMsg(paste0("  Initiation points: ", initPointsFile, "\n"))
  logMsg(paste0("  Non-initiation points ratio: ", noninitRatio, "\n"))
  logMsg(paste0("  Buffer radius: ", bufferRadius, "\n"))
  logMsg(paste0("  Buffer extraction method: ", bufferExtractionMethod, "\n"))
  logMsg(paste0("  Initiation range expansion: ", initRangeExpansion, "%\n"))
  logMsg(paste0("  Iterations: ", iterationsCount, "\n"))
  logMsg(paste0("  Testing proportion: ", testingProportion, "%\n"))
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
  initRange <- createInitiationRange(
    varsRaster,
    initBuffers,
    initRangeExpansion
  )
  
  # Log initiation range matrix
  logMsg("INITIATION RANGES:\n")
  logObj(initRange)
  logMsg("\n")
  
  # Identify cells in the study region that have variable values within their 
  # initiation ranges
  analysisRegionMask <- createAnalysisRegionMask(
    varsRaster,
    initRange
  )
  
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
  
  # Generate non-initiation buffers
  noninitBuffers <- generateNoninitiationBuffers(
    noninitBuffersCount,
    noninitRegion,
    bufferRadius
  )
  
  # Create non-initiation buffer sets ------------------------------------------
  
  # Create a static set of training and testing non-initiation buffers to use every iteration
  testingNoninitCount <- floor(length(noninitBuffers) * (testingProportion / 100))
  testingNoninitIndices <- sample(seq_along(noninitBuffers), size = testingNoninitCount)
  trainingNoninitIndices <- setdiff(seq_along(noninitBuffers), testingNoninitIndices)
  
  trainingNoninitBuffers <- noninitBuffers[trainingNoninitIndices]
  testingNoninitBuffers <- noninitBuffers[testingNoninitIndices]
  
  # Perform cross-validation ---------------------------------------------------
  
  # Place to store iteration model AUC values
  iterationsAucValues <- c()
  
  # Place to store iteration model error rates
  iterationsErrorRates <- data.frame(
    rep(NA, iterationsCount),
    rep(NA, iterationsCount),
    rep(NA, iterationsCount)
  )
  names(iterationsErrorRates) <- c("OOB", "initiation", "non-initiation")
  
  # Calculate how many initiation buffers should be used for testing per iteration
  testingInitBuffersCount <- floor(length(initBuffers) * (testingProportion / 100))
  
  for (i in seq_len(iterationsCount)) {
    ## Create initiation buffer sets -------------------------------------------
    
    # Create testing initiation set
    testingInitIndices <- sample(seq_along(initBuffers), size = testingInitBuffersCount)
    testingInitBuffers <- initBuffers[testingInitIndices]
    
    # Create training initiation set
    trainingInitIndices <- setdiff(seq_along(initBuffers), testingInitIndices)
    trainingInitBuffers <- initBuffers[trainingInitIndices]
    
    ## Create model ------------------------------------------------------------
    
    # Create training dataset
    trainingData <- extractBufferValues(
      varsRaster,
      trainingInitBuffers,
      trainingNoninitBuffers,
      bufferExtractionMethod
    )
    
    coordsCols <- names(trainingData) %in% c("x", "y")  
    trainingData <- trainingData[,!coordsCols]
    
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
      testingInitBuffers,
      testingNoninitBuffers,
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
    
    # Log AUC value
    auc <- rocStats$auc@y.values[[1]]
    logMsg(paste0("TESTING AUC: ", round(auc, digits = 7), "\n"))
    logMsg("\n")
    
    # Record iteration statistics
    iterationsAucValues <- c(iterationsAucValues, auc)
    iterationsErrorRates[i,] <- as.data.frame(rfModel$err.rate)[rfModel$ntree,]
  }
  
  # Summarize iterations -------------------------------------------------------
  
  logMsg("SUMMARY --------------------------------------------\n\n")
  
  # Summarize AUC values
  aucValuesSummary <- summarizeVector(iterationsAucValues)
  rownames(aucValuesSummary) <- "AUC"
  
  # Summarize error rates
  errorRatesSummary <- summarizeDataFrame(iterationsErrorRates)
  
  # Log summary of AUC values
  logMsg("TESTING AUC:\n")
  logObj(aucValuesSummary)
  logMsg("\n")
  
  # Log summary of error rates
  logMsg("TESTING ERROR RATES:\n")
  logObj(errorRatesSummary)
  
}

# Entrypoint for ArcGIS script tool
tool_exec <- function(in_params, out_params) {
  
  assessInitiationPointSusceptibility(
    refRasterFile          = in_params[[1]],
    varRasterFiles         = in_params[[2]],
    initPointsFile         = in_params[[3]],
    noninitRatio           = in_params[[4]],
    bufferRadius           = in_params[[5]],
    bufferExtractionMethod = in_params[[6]],
    initRangeExpansion     = in_params[[7]],
    iterationsCount        = in_params[[8]],
    testingProportion      = in_params[[9]],
    outputDir              = in_params[[10]]
  )
  
  return(out_params)
  
}
 