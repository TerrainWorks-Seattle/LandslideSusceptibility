#' @title Create non-initiation raster
#'
#' @description Creates a raster that defines where non-initiation points can be
#' generated. This includes cells not already in initiation sites and whose 
#' variable values are all within their respective initiation ranges.
#'
#' @param raster      A SpatRaster of explanatory variables
#' @param initRange   A matrix that holds the min & max initiation limits 
#'                    of each explanatory variable
#' @param initBuffers A SpatVector of initiation buffers
#' 
#' @return A SpatRaster where any non-NA cell is a viable location for a 
#' non-initiation point.
#' @export

createNoninitiationRaster <- function(raster, initRange, initBuffers) {
  
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
  initRaster <- terra::app(initRaster, fun = "prod")
  
  # NA-out cells within initiation buffers
  initCells <- terra::extract(initRaster, initBuffers, cells = TRUE)$cell
  initRaster[initCells] <- NA
  
  return(initRaster)
  
}