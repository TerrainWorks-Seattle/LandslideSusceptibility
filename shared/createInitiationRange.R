# Creates a matrix that holds the min/max initiation limits of each 
# explanatory variable.

createInitiationRange <- function(
  varsRaster,              # A SpatRaster of explanatory variables
  initiationBuffers,       # A SpatVector of initiation buffers
  initiationRangeExpansion # Percent to reduce each range min and increase each range max by 
) {
  # Extract all variable values from initiation buffers
  initiationValues <- terra::extract(varsRaster, initiationBuffers)
  
  # For each initiation buffer, find each variable's maximum value
  regionMaxVarValues <- aggregate(. ~ ID, data = initiationValues, max)
  regionMaxVarValues$ID <- NULL # Remove "ID" column
  initiationValues$ID   <- NULL # Remove "ID" column
  
  # Find the min and max of the maximum variable values in each buffer 
  initiationMinValues <- lapply(regionMaxVarValues, min)
  initiationMaxValues <- lapply(regionMaxVarValues, max)
  
  # Slightly expand the initiation range
  expansionRatio <- initiationRangeExpansion / 100  
  initiationMinValues <- lapply(initiationMinValues, function(x) x - expansionRatio * x)
  initiationMaxValues <- lapply(initiationMaxValues, function(x) x + expansionRatio * x)
  
  # Create a matrix that holds each variable's initiation range
  initiationRange <- matrix(rep(NA, 4), nrow = 2)
  colnames(initiationRange) <- c("min", "max")
  rownames(initiationRange) <- names(varsRaster)
  for (varName in names(varsRaster)) {
    initiationRange[varName, "min"] <- initiationMinValues[[varName]]
    initiationRange[varName, "max"] <- initiationMaxValues[[varName]]
  }
  
  return(initiationRange)
}