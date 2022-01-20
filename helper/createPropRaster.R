createPropRaster <- function(probRaster) {
  
  # Order probability values from least to greatest 
  probValues <- terra::values(probRaster)
  valueOrder <- order(probValues)
  
  # Calculate cumulative sum of all non-NA probability values
  naIndices <- which(is.na(probValues[valueOrder]))
  probCumSum <- if (length(naIndices) == 0) {
    cumsum(probValues[valueOrder])
  } else {
    cumsum(probValues[valueOrder[1:naIndices[1]-1]])
  }
  
  # Calculate proportion values
  totalProbSum <- probCumSum[length(probCumSum)]
  propValues <- probCumSum / totalProbSum
  
  # Add back NA values if any were present
  propValues <- c(propValues, rep(NA, length(naIndices)))
  
  # Re-order proportion values to align back with their corresponding cells
  propValues <- propValues[order(valueOrder)]
  
  # Create proportion raster
  propRaster <- terra::rast(probRaster)
  terra::values(propRaster) <- propValues
  
  propRaster
}
