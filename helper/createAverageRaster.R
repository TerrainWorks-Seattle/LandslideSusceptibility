#' @title Create average raster
#' 
#' @description Creates an average raster from a set of raster files. Good for 
#' averaging many large raster files that might not fit into memory all at once.
#' 
#' @param rasterFiles a vector of raster files to average together
#' 
#' @returns a raster object.
#' 
#' @example
#' \donttest{
#' }

createAverageRaster <- function(rasterFiles) {
  
  if (length(rasterFiles) == 0)
    return(NULL)
  
  # Average all rasters
  sumRaster <- terra::rast(rasterFiles[1])
  if (length(rasterFiles) > 1) {
    for (i in 2:length(rasterFiles)) {
      sumRaster <- sumRaster + terra::rast(rasterFiles[i])
    }
  }
  avgRaster <- sumRaster / length(rasterFiles)
  
  return(avgRaster)
  
}
