#' @title Create a proportion raster
#' 
#' @description Creates a proportion raster from a probability raster.
#' 
#' @param probRaster A raster containing probability values from 0.0 to 1.0
#' 
#' @examples
#' \donttest{
#' r1 <- terra::rast(nrow = 3, ncol = 4)
#' terra::values(r1) <- runif(terra::ncell(r1))
#' terra::values(r1)[9] <- NA
#' terra::values(r1)[7] <- NA
#' p1 <- createProportionRaster(r1)
#' 
#' r2 <- terra::rast("E:/NetmapData/Scottsburg/prob_avg.tif")
#' p2 <- createProportionRaster(r2)
#' }

createProportionRaster <- function(probRaster) {
  
  # Calculate the cumulative sum for each cell in order from least to greatest 
  # probability
  probValues <- terra::values(probRaster)
  valueOrder <- order(probValues)
  naIndices <- which(is.na(probValues[valueOrder]))
  probCumSum <- if (length(naIndices) == 0) {
    cumsum(probValues[valueOrder])
  } else {
    # Avoid summing NA values if present
    cumsum(probValues[valueOrder[1:naIndices[1]-1]])
  }
  
  # Calculate proportions by dividing cumsums by the total probability sum
  totalProbSum <- probCumSum[length(probCumSum)]
  propValues <- probCumSum / totalProbSum
  
  # Add back NA values if any were present
  propValues <- c(propValues, rep(NA, length(naIndices)))
  
  # Re-order proportion values to align with their corresponding cells
  propValues <- propValues[order(valueOrder)]
  
  # Create a proportion raster
  propRaster <- terra::rast(probRaster)
  terra::values(propRaster) <- propValues

  return(propRaster)
  
}
