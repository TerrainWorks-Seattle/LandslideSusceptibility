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