refRasterFile          <- "E:/NetmapData/Scottsburg/elev_scottsburg.flt"
varRasterFiles         <- list(
                            "E:/NetmapData/Scottsburg/grad_30.tif",
                            "E:/NetmapData/Scottsburg/plan_30.tif"
                          )
initPointsFile         <- "E:/NetmapData/Scottsburg/Scottsburg_Upslope.shp"
noninitRatio           <- 1.5
bufferRadius           <- 20
bufferExtractionMethod <- "center cell"
initRangeExpansion     <- 0
noninitSets            <- 10 # How many sets of non-init points to generate and CV with
cvRepetitions          <- 2  # How many times should k-fold CV be repeated?
testingProportion      <- 10
outputDir              <- "E:/NetmapData/Scottsburg"

# Package dependencies
# mlr3
# mlr3learners
# mlr3spatiotempcv
# mlr3viz
# precrec
# ranger

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
if (!(bufferExtractionMethod %in% c("all cells", "center cell", "max gradient 
      cell", "max plan cell")))
  stop("Buffer extraction method must be one of 'all cells', 'center cell', 
       'max gradient cell', or 'max plan cell'.")

# Validate number of iterations
if (iterations < 1)
  stop("Iterations cannot be fewer than 1.")

# Validate testing proportion
if (testingProportion <= 0 || testingProportion >= 100)
  stop("Testing proportion must be between 0% and 100%.")

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

# Perform cross-validation -----------------------------------------------------

for (i in seq_len(noninitSets)) {
  
  # Generate non-initiation buffers --------------------------------------------
  
  noninitBuffers <- generateNoninitiationBuffers(
    noninitBuffersCount,
    noninitRegion,
    bufferRadius
  )
  
  # Create full landslide dataset ----------------------------------------------
  
  landslideData <- extractBufferValues(
    varsRaster,
    initBuffers,
    noninitBuffers,
    bufferExtractionMethod
  )
  
  # Define cross-validation method ---------------------------------------------
  
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
  
  # Create a version of the variables raster that only keeps cell values within 
  # the analysis area
  analysisVarsRaster <- terra::mask(varsRaster, analysisRegionMask)
  
  # Set up a resampling method for repeated k-fold spatial cross-validation
  resampling <- mlr3::rsmp("repeated_spcv_coords", folds = 5, repeats = cvRepetitions)
  
  # Perform the resampling method on the task
  resampling <- resampling$instantiate(task)
  
  # --------------------------------------------------------------------
  
  for (j in seq_len(resampling$iters)) {
    
    ## Train model -------------------------------------------------------------
    
    # Get training data
    trainingIndices <- resampling$train_set(j)
    trainingData <- landslideData[trainingIndices]
    
    # Remove columns for coordinates
    coordsCols <- names(trainingData) %in% c("x", "y")  
    trainingData <- trainingData[,!coordsCols]
    
    # Train random forest model
    rfModel <- randomForest::randomForest(
      formula = class ~ .,
      data = trainingData
    )
    
    ## Evaluate model ----------------------------------------------------------
    
    # Get testing data
    testingIndices <- resampling$test_set(j)
    testingData <- landslideData[testingIndices]
    
    # Remove columns for coordinates
    coordsCols <- names(testingData) %in% c("x", "y")  
    testingData <- testingData[,!coordsCols]
    
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
    
  }
  
}

# # Set up a random forest learner
# learner <- mlr3::lrn("classif.ranger", predict_type = "prob")
# 
# # Plot each fold's train/test points distribution for iteration 1
# #mlr3spatiotempcv::autoplot(resampling, task, fold_id = c(1:5), repeats_id = 1)
# 
# # Perform resampling method with the random forest learner and the landslides task
# result <- mlr3::resample(
#   task = task,
#   learner = learner,
#   resampling = resampling
# )
#
# # Inspect results --------------------------------------------------------
# 
# # Model 1 ROC plot
# mlr3viz::autoplot(result$predictions()[[4]], type = "roc")
# 
# # Model 1 predictions
# result$predictions()[[1]]
# 
# # Model 1 confusion matrix
# result$predictions()[[1]]$confusion
# 
# # Calculate average classification error
# result$aggregate(measures = mlr3::msr("classif.ce"))
# 
# # Calculate average area under ROC curve
# result$aggregate(measures = mlr3::msr("classif.auc"))
# 
# # Calculate average accuracy (correct predictions out of all predictions)
# result$aggregate(measures = mlr3::msr("classif.acc"))
