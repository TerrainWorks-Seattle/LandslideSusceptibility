#' @title Create analysis region mask
#'
#' @description Creates a raster mask that identifies area(s) with landslide 
#' initiation characteristics. Cells with variable values outside their 
#' respective initiation ranges are marked with NA.
#'
#' @param raster    A SpatRaster of explanatory variables
#' @param initRange A matrix that holds the min & max initiation limits of each 
#'                  explanatory variable
#' 
#' @return A SpatRaster with NA cells anywhere outside the analysis region.

createAnalysisRegionMask <- function(raster, initRange) {
  
  initRaster <- terra::copy(raster)
  for (varName in names(initRaster)) {
    varRaster <- initRaster[[varName]]
    
    # Get variable value limits
    minInitValue <- initRange[varName, "min"]
    maxInitValue <- initRange[varName, "max"]
    
    # NA-out cells with values outside variable initiation range
    varInitRaster <- terra::app(varRaster, function(x) {
      ifelse(x < minInitValue | x > maxInitValue, NA, x)
    })
    
    # Update the raster in the input raster stack
    initRaster[[varName]] <- varInitRaster
  }
  
  # NA-out cells with ANY variable value outside its initiation range
  analysisRegionMask <- terra::app(initRaster, fun = "all")
  
  return(analysisRegionMask)
  
}