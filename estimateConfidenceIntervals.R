
estimateConfidenceIntervals <- function(
  refRasterFile,
  varRasterFiles,
  initPointsFile,
  noninitRatio,
  bufferRadius,
  bufferExtractionMethod,
  initRangeExpansion,
  folds,
  repetitions,
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
  
  source("./shared/createInitiationRange.R")
  source("./shared/createAnalysisRegionMask.R")
  source("./shared/generateNoninitiationBuffers.R")
  source("./shared/extractBufferValues.R")
  source("./shared/createProportionRaster.R")
  
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
  
  # Validate number of repetitions
  if (repetitions < 1)
    stop("Repetitions cannot be fewer than 1.")
  
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
  
  # Plot variable distributions ------------------------------------------------
  
  # Get requested values from all initiation and non-initiation buffers
  landslideData <- extractBufferValues(
    varsRaster,
    initBuffers,
    noninitBuffers,
    bufferExtractionMethod
  )
  
  # Create cross-validation sets ----------------------------------------------
  
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
    folds = folds,
    repeats = repetitions
  )
  
  totalIterations <- resampling$iters
  
  # Generate each CV iteration's training/testing sets using the resampling method
  cvSets <- resampling$instantiate(task)
  
  # Perform cross-validation ---------------------------------------------------
  
  probRasterList <- list()
  
  # Create a version of the variables raster that only keeps the cells within 
  # the analysis region
  analysisVarsRaster <- terra::mask(varsRaster, analysisRegionMask)
  
  for (i in seq_len(totalIterations)) {
    ## Get iteration training/testing data ------------------------------------
    
    trainingSet <- cvSets$train_set(i)
    testingSet <- cvSets$test_set(i)
    
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
  probConfidenceRaster <- createConfidenceIntervalRaster(probsRaster, p = 0.05)
  propConfidenceRaster <- createConfidenceIntervalRaster(propsRaster, p = 0.05)
  
  terra::writeRaster(probConfidenceRaster, paste0(outputDir, "/probconf.tif"))
  terra::writeRaster(propConfidenceRaster, paste0(outputDir, "/propconf.tif"))
}

createConfidenceIntervalRaster <- function(x, p = 0.05) {
  n <- terra::nlyr(x)
  degreesFreedom <- n - 1
  tScore <- qt(p = p, df = degreesFreedom, lower.tail = FALSE)
  
  sd <- terra::app(x, fun = "sd") # Standard deviation
  se <- sd / sqrt(n)              # Standard error
  me <- se * tScore               # Margin of error
  
  mean <- terra::mean(x)          # Mean
  lower <- mean - me              # Interval lower bound
  upper <- mean + me              # Interval upper bound
  
  interval <- c(lower, upper)
  names(interval) <- c("lower", "upper")
  
  return(interval)
}

if (FALSE) {
  estimateConfidenceIntervals(
    refRasterFile          = "E:/NetmapData/Scottsburg/elev_scottsburg.flt",
    varRasterFiles         = list(
                               "E:/NetmapData/Scottsburg/grad_30.tif",
                               "E:/NetmapData/Scottsburg/plan_30.tif",
                               "E:/NetmapData/Scottsburg/pca_scott.flt"
                             ),   
    initPointsFile         = "E:/Netmapdata/Scottsburg/Scottsburg_Upslope.shp",
    noninitRatio           = 1.5,
    bufferRadius           = 20,
    bufferExtractionMethod = "center cell",
    initRangeExpansion     = 0,
    folds                  = 5,
    repetitions            = 4,
    outputDir              = "E:/NetmapData"
  )
}
