#' @title Generate non-initiation buffers
#' 
#' @description Generates a SpatVector of non-initiation buffers in a given 
#' area.
#' 
#' @param count  The number of non-initiation buffers to generate
#' @param region A SpatRaster with NA cells anywhere a non-initiation buffer
#'               cannot be located
#' @param radius The radius of each non-initiation buffer
#' 
#' @return A SpatVector of non-initiation buffers.
#' @export

generateNoninitiationBuffers <- function(count, region, radius) {
  
  # NOTE: terra::spatSample() sometimes generates less than the requested 
  # number of points if the sample raster has a lot of NAs. This is remedied 
  # by repeatedly requesting a larger and larger number of points until
  # enough have been generated, then subsetting those for the correct amount.
  
  currentRequest <- count
  hasGeneratedEnough <- FALSE
  
  while (!hasGeneratedEnough) {
    # Sample points anywhere that fits initiation conditions but recorded no landslides
    noninitPoints <- terra::spatSample(
      region,
      size = currentRequest,
      na.rm = TRUE,
      as.points = TRUE,
      warn = FALSE
    )
    
    # Test if enough non-initiation points have been generated
    if (length(noninitPoints) >= count) {
      # If so, subset and exit loop
      noninitPoints <- noninitPoints[seq_len(count)]
      hasGeneratedEnough <- TRUE
    } else {
      # If not, double the next request
      currentRequest <- currentRequest * 2
    }
  }
  
  # Create a buffer around each non-initiation point
  noninitBuffers <- if (radius > 0) {
    terra::buffer(noninitPoints, width = radius)
  } else {
    noninitPoints
  }
  
  return(noninitBuffers)
  
}