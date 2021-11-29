#' @title Estimate confidence intervals
#' 
#' @description Generates a confidence interval raster from a series of 
#' landslide initiation probability rasters. These are produced by models 
#' created
#' 
#' @param refRasterFile          A raster file to use as a grid reference.
#' @param varRasterFiles         A list of raster files to use as explanatory 
#'                               variables.
#' @param initPointsFile         A shapefile of initiation points.
#' @param noninitProportion      The proportion of non-initiation sites to 
#'                               initiation sites.
#' @param bufferRadius           The radius of site buffers.
#' @param bufferExtractionMethod The method used to select values from site 
#'                               buffers for training/testing. Either: 
#'                               "all cells", "center cell", 
#'                               "max gradient cell", or "max plan cell".
#' @param initRangeExpansion     The proportion (in %) to expand the initiation 
#'                               range of each variable by.
#' @param repetitionsCount       The number of times to repeat k-fold cross-
#'                               validation.
#' @param foldsCount             The number of folds to use for k-fold cross-
#'                               validation.
#' @param pValue                 P-value of confidence intervals.
#' @param outputDir              The directory to write output files to.
#' 
#' @example
#' \dontest{
#' estimateConfidenceIntervals(
#'   refRasterFile          = "E:/NetmapData/Scottsburg/elev_scottsburg.flt",
#'   varRasterFiles         = list(
#'     "E:/NetmapData/Scottsburg/grad_30.tif",
#'     "E:/NetmapData/Scottsburg/plan_30.tif"
#'   ),   
#'   initPointsFile         = "E:/Netmapdata/Scottsburg/Scottsburg_Upslope.shp",
#'   noninitProportion      = 1,
#'   bufferRadius           = 20,
#'   bufferExtractionMethod = "center cell",
#'   initRangeExpansion     = 0,
#'   repetitionsCount       = 1,
#'   foldsCount             = 5,
#'   pValue                 = 0.05,
#'   outputDir              = "E:/NetmapData/Scottsburg"
#' )
#' }

estimateConfidenceIntervals <- function(
  refRasterFile,
  varRasterFiles,
  initPointsFile,
  noninitProportion,
  bufferRadius,
  bufferExtractionMethod,
  initRangeExpansion,
  repetitionsCount,
  foldsCount,
  pValue,
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
  source("./helper/createProportionRaster.R")
  source("./helper/createConfidenceIntervalRaster.R")
  
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
  
  # Validate non-initiation proportion
  if (noninitProportion <= 0)
    stop("Non-initiation points proportion must be greater than 0.")
  
  # Validate buffer radius
  if (bufferRadius < 0)
    stop("Buffer radius must be greater than or equal to 0.")
  
  # Validate buffer extraction method
  if (!(bufferExtractionMethod %in% c("all cells", "center cell", 
                                      "max gradient cell", "max plan cell")))
    stop("Buffer extraction method must be one of 'all cells', 'center cell', 
       'max gradient cell', or 'max plan cell'.")
  
  # Validate number of repetitions
  if (repetitionsCount < 1)
    stop("Repetitions must be greater than or equal to 1.")
  
  # Validate number of folds
  if (foldsCount < 2)
    stop("Folds must be greater than 1.")
  
  # Validate output directory
  if (!file.exists(outputDir))
    stop("Output directory not found: '", outputDir, "'.")
  
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
  
  # Identify cells in the study region that have variable values within their 
  # initiation ranges
  analysisRegionMask <- createAnalysisRegionMask(
    varsRaster,
    initRange
  )
  
  # Define the region where non-initiation points can be generated
  
  # NOTE: initiation buffers and non-initiation buffers must not overlap. This 
  # can be avoided by making sure each non-initiation point is generated at 
  # least 2 buffer-radius-lengths away from any initiation point 
  
  # Double the size of the initiation buffers
  expInitBuffers <- terra::buffer(initPoints, width = bufferRadius * 2)
  
  # Remove expanded initiation buffers from the viable non-initiation region
  noninitRegion <- terra::copy(analysisRegionMask)
  initCellIndices <- terra::extract(noninitRegion, expInitBuffers, cells = TRUE)$cell
  noninitRegion[initCellIndices] <- NA
  
  # Determine how many non-initiation buffers to generate
  noninitBuffersCount <- ceiling(length(initPoints) * noninitProportion)
  
  # Generate non-initiation buffers
  noninitBuffers <- generateNoninitiationBuffers(
    noninitBuffersCount,
    noninitRegion,
    bufferRadius
  )
  
  # Create cross-validation sets -----------------------------------------------
  
  # Get requested values from all initiation and non-initiation buffers
  landslideData <- extractBufferValues(
    varsRaster,
    initBuffers,
    noninitBuffers,
    bufferExtractionMethod
  )
  
  # Set up a machine learning task for the landslide data
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
  
  # Set up a resampling method for repeated k-fold spatial cross-validation
  resampling <- mlr3::rsmp(
    "repeated_spcv_coords",
    folds = foldsCount,
    repeats = repetitionsCount
  )
  
  totalIterations <- resampling$iters
  
  # Generate each CV iteration's training/testing sets using the resampling method
  resampling <- resampling$instantiate(task)
  
  # Perform cross-validation ---------------------------------------------------
  
  probRasterList <- list()
  
  # Create a version of the variables raster that only keeps the cells within 
  # the analysis region
  analysisVarsRaster <- terra::mask(varsRaster, analysisRegionMask)
  
  # Remove coordinates columns from landslide dataset
  coordsCols <- names(landslideData) %in% c("x", "y")  
  landslideData <- landslideData[,!coordsCols]
  
  for (i in seq_len(totalIterations)) {
    
    ## Get iteration training/testing data ------------------------------------
    
    trainingSet <- resampling$train_set(i)
    testingSet <- resampling$test_set(i)
    
    # Remove coordinates from dataset
    coordsCols <- names(landslideData) %in% c("x", "y")  
    trainingData <- landslideData[trainingSet, !coordsCols]
    testingData <- landslideData[testingSet, !coordsCols]
    
    ## Create model ------------------------------------------------------------
    
    # Train a new random forest model
    rfModel <- randomForest::randomForest(
      formula = class ~ .,
      data = trainingData
    )
    
    ## Generate probability raster ---------------------------------------------
    
    probRaster <- terra::predict(
      analysisVarsRaster,
      rfModel,
      na.rm = TRUE,
      type = "prob"
    )[["initiation"]]
    
    probRasterList[[i]] <- probRaster
    
  }
  
  # Create confidence interval rasters -----------------------------------------

  # Create landslide proportion rasters
  propRasterList <- lapply(probRasterList, function(x) createProportionRaster(x))
  
  # Combine individual rasters into single multi-layered rasters
  probsRaster <- terra::rast(probRasterList)
  propsRaster <- terra::rast(propRasterList)
  
  # Create confidence interval rasters
  probConfRaster <- createConfidenceIntervalRaster(probsRaster, p = pValue)
  propConfRaster <- createConfidenceIntervalRaster(propsRaster, p = pValue)
  
  # Save confidence interval rasters
  terra::writeRaster(probConfRaster, paste0(outputDir, "/prob_conf.tif"), overwrite = TRUE)
  terra::writeRaster(propConfRaster, paste0(outputDir, "/prop_conf.tif"), overwrite = TRUE)
}
