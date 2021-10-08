# Creates a raster with NA cells anywhere the input variable rasters don't
# fall within their respective initiation ranges or are inside initiation
# sites.

createNoninitiationRaster <- function(
  varsRaster,       # A SpatRaster of explanatory variables
  initiationRange,  # A matrix that holds the min/max initiation limits of each explanatory variable 
  initiationBuffers # A SpatVector of initiation buffers
) {
  initiationRaster <- terra::copy(varsRaster)
  for (varName in names(initiationRaster)) {
    varRaster <- initiationRaster[[varName]]
    
    # Get variable value limits
    minInitiationValue <- initiationRange[varName, "min"]
    maxInitiationValue <- initiationRange[varName, "max"]
    
    # NA-out cells with values outside variable initiation range
    varInitiationRaster <- terra::app(varRaster, function(x) {
      ifelse(x < minInitiationValue | x > maxInitiationValue, NA, x)
    })
    
    # Update the raster in the input raster stack
    initiationRaster[[varName]] <- varInitiationRaster
  }
  
  # NA-out cells with variable values outside their initiation range
  initiationRaster <- terra::app(initiationRaster, fun = "prod")
  
  # NA-out cells within initiation buffers
  initiationCells <- terra::extract(initiationRaster, initiationBuffers, cells = TRUE)$cell
  initiationRaster[initiationCells] <- NA
  
  return(initiationRaster)
}