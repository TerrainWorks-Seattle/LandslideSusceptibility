#' @title Create initiation range
#'
#' @description Finds the initiation range for each raster variable.
#'
#' @param raster    A SpatRaster of explanatory variables
#' @param initBuffers   A SpatVector of initiation buffers
#' @param expansion Percent to reduce each range min and increase each range 
#' max by 
#'
#' @return A matrix that holds the min & max initiation limits of each raster 
#' variable.

createInitiationRange <- function(raster, initBuffers, expansion) {
  
  # Extract all variable values from initiation buffers
  initiationValues <- terra::extract(raster, initBuffers)
  
  # For each initiation buffer, find each variable's maximum value
  regionMaxVarValues <- aggregate(. ~ ID, data = initiationValues, max)
  regionMaxVarValues$ID <- NULL # Remove "ID" column
  initiationValues$ID   <- NULL # Remove "ID" column
  
  # Find the min and max of the maximum variable values in each buffer 
  initiationMinValues <- lapply(regionMaxVarValues, min)
  initiationMaxValues <- lapply(regionMaxVarValues, max)
  
  # Expand each initiation range
  initiationMinValues <- lapply(initiationMinValues, function(x) x - (expansion / 100) * x)
  initiationMaxValues <- lapply(initiationMaxValues, function(x) x + (expansion / 100) * x)
  
  # Create a matrix that holds each variable's initiation range
  varCount <- terra::nlyr(raster)
  initiationRange <- matrix(rep(NA, varCount * 2), nrow = varCount)
  colnames(initiationRange) <- c("min", "max")
  rownames(initiationRange) <- names(raster)
  for (varName in names(raster)) {
    initiationRange[varName, "min"] <- initiationMinValues[[varName]]
    initiationRange[varName, "max"] <- initiationMaxValues[[varName]]
  }
  
  return(initiationRange)
  
}