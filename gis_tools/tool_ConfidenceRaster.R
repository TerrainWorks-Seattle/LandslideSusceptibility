#' @title Create a confidence raster
#' 
#' @description Creates a confidence raster by calculating the confidence 
#' interval for each cell across a list of input rasters.
#' 
#' @param rasterFiles A list of raster files.
#' @param p           The p-value of each cell's confidence interval.
#' @param outputFile  The raster file to output.

createConfidenceRaster <- function(
  rasterFiles,
  p,
  outputFile
) {

  
  # Validate parameters --------------------------------------------------------
  
  # Validate raster files
  lapply(rasterFiles, function(file) {
    if (!file.exists(file))
      stop(paste0("Raster file not found: '", file, "'."))
  })
  
  # Validate p-value
  if (p <= 0 || p >= 1)
    stop("P-value must be 0 < p < 1.")

  # Create confidence interval raster ------------------------------------------
  
  # Load rasters
  rasterList <- lapply(rasterFiles, function(file) terra::rast(file))
  
  # Combine individual rasters into a single multi-layered raster
  raster <- terra::rast(rasterList)
  
  n <- terra::nlyr(raster)
  degreesFreedom <- n - 1
  tScore <- qt(p = p, df = degreesFreedom, lower.tail = FALSE)
  
  sd <- terra::app(raster, fun = "sd") # Standard deviation
  se <- sd / sqrt(n)                   # Standard error
  me <- se * tScore                    # Margin of error
  
  mean <- terra::mean(raster)          # Mean
  lower <- mean - me                   # Interval lower bound
  upper <- mean + me                   # Interval upper bound
  
  # Combine upper and lower bound layers into a single raster
  confRaster <- c(lower, upper)
  names(confRaster) <- c("lower", "upper")
  
  # Save confidence raster
  terra::writeRaster(confRaster, outputFile, overwrite = TRUE)
  
}

# Entrypoint for ArcGIS script tool
tool_exec <- function(in_params, out_params) {
  
  createConfidenceRaster(
    rasterFiles = in_params[[1]],
    p           = in_params[[2]],
    outputFile  = out_params[[1]]
  )
  
  return(out_params)
  
}
