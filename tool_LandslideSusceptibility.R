#' @title Predict landslide susceptibility
#'
#' @description Creates a landslide susceptibility raster using a random forest
#' model trained on initiation sites and topographic variable rasters. The value
#' of each cell represents its 0.0-1.0 prediction score of matching landslide
#' initiation conditions.
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
#' @param iterations             How many models to create (trained on different
#'                               sets of non-initiation sites)
#' @param outputDir              The directory to write output files to
#' 
#' @example
#' \donttest{
#' predictLandslideSusceptibility(
#'   refRasterFile = "E:/NetmapData/Scottsburg/elev_scottsburg.flt",
#'   varRasterFiles <- list(
#'     "E:/NetmapData/Scottsburg/grad_30.tif",
#'     "E:/NetmapData/Scottsburg/plan_30.tif",
#'     "E:/Netmapdata/Scottsburg/pca_scott.flt"
#'   ),
#'   initPointsFile = "E:/NetmapData/Scottsburg/Scottsburg_Upslope.shp",
#'   noninitRatio = 1.5,
#'   bufferRadius = 20,
#'   bufferExtractionMethod = "center cell",
#'   initRangeExpansion = 50,
#'   iterations = 1,
#'   outputDir = "E:/NetmapData/Scottsburg"
#' )
#' }

predictLandslideSusceptibility <- function(
  refRasterFile,
  varRasterFiles,
  initPointsFile,
  noninitRatio,
  bufferRadius,
  bufferExtractionMethod,
  initRangeExpansion,
  iterations,
  outputDir
) {
  
  # Load helper functions ------------------------------------------------------
  
  # NOTE: running an R ArcGIS script tool automatically sets the current working
  # directory to the folder that holds the tool's script file, so any sourced 
  # file paths must be relative to that folder.
  
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

  # Validate number of iterations
  if (iterations < 1)
    stop("Iterations cannot be fewer than 1.")
  
  # Validate output directory
  if (!file.exists(outputDir))
    stop("Output directory not found: '", outputDir, "'.")
  
  # Set up logging -------------------------------------------------------------
  
  logFilename <- "landslide_susceptibility.txt"
  file.create(paste0(outputDir, "/", logFilename))
  
  logObj <- function(obj) {
    capture.output(obj, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  logMsg <- function(msg) {
    cat(msg, file = paste0(outputDir, "/", logFilename), append = TRUE)
  }
  
  # Log parameter values
  logMsg(paste0(getwd(), "\n"))
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
  logMsg(paste0("  Iterations: ", iterations, "\n"))
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
  
  # Prepare to generate non-initiation buffers ---------------------------------
  
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

  # Determine how many to generate
  noninitBuffersCount <- ceiling(length(initPoints) * noninitRatio)
    
  # Generate a set of initiation probability rasters ---------------------------
  
  # Place to store initation probability rasters
  initProbRasterList <- list()
  
  # Create a version of the variables raster that only keeps cell values within 
  # the analysis area
  analysisAreaVarsRaster <- terra::mask(varsRaster, analysisRegionMask)
  
  for (i in seq_len(iterations)) {
    ## Generate non-initiation buffers -----------------------------------------
    
    # Generate non-initiation buffers
    noninitBuffers <- generateNoninitiationBuffers(
      noninitBuffersCount,
      noninitRegion,
      bufferRadius
    )
    
    ## Create model ------------------------------------------------------------
    
    # Create training dataset
    trainingData <- extractBufferValues(
      varsRaster,
      initBuffers,
      noninitBuffers,
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

    initProbRaster <- terra::predict(
      analysisAreaVarsRaster,
      rfModel,
      na.rm = TRUE,
      type = "prob"
    )[["initiation"]]

    initProbRasterList[[i]] <- initProbRaster
  }
  
  # Summarize initiation probability rasters -----------------------------------
  
  if (length(initProbRasterList) > 1) {
    # Combine probability rasters into a single multi-layer raster
    initProbRaster <- terra::rast(initProbRasterList)

    # Generate summary initiation probability rasters
    avgInitProbRaster <- terra::app(initProbRaster, fun = "mean")
    minInitProbRaster <- terra::app(initProbRaster, fun = "min")
    maxInitProbRaster <- terra::app(initProbRaster, fun = "max")

    # Save summary rasters
    terra::writeRaster(avgInitProbRaster, paste0(outputDir, "/prob_avg.tif"), overwrite = TRUE)
    terra::writeRaster(minInitProbRaster, paste0(outputDir, "/prob_min.tif"), overwrite = TRUE)
    terra::writeRaster(maxInitProbRaster, paste0(outputDir, "/prob_max.tif"), overwrite = TRUE)
  } else {
    # Save single probability raster
    terra::writeRaster(initProbRasterList[[1]], paste0(outputDir, "/prob.tif"), overwrite = TRUE)
  }
  
}

# ArcGIS script tool entrypoint
tool_exec <- function(in_params, out_params) {
  
  predictLandslideSusceptibility(
    refRasterFile          = in_params[[1]],
    varRasterFiles         = in_params[[2]],
    initPointsFile         = in_params[[3]],
    noninitRatio           = in_params[[4]],
    bufferRadius           = in_params[[5]],
    bufferExtractionMethod = in_params[[6]],
    initRangeExpansion     = in_params[[7]],
    iterations             = in_params[[8]],
    outputDir              = in_params[[9]]
  )
  
  return(out_params)
  
}