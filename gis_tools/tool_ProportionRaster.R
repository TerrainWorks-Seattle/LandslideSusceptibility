#' @title Create a proportion raster
#' 
#' @description Creates a proportion raster from a probability raster, allowing 
#' you to see where a given proportion of landslides will likely occur.
#' 
#' @param probRasterFile A raster containing probability values from 0.0 to 1.0.
#' @param outputFile     The raster file to output.
#' 
#' @examples
#' \donttest{
#' createProportionRaster(
#'   "E:/NetmapData/Scottsburg/avg_prob.tif",
#'   "E:/NetmapData/Scottsburg/proportion.tif"
#' )
#' }

createProportionRaster <- function(
  probRasterFile,
  outputFile
) {
  
  # Install dependencies -------------------------------------------------------
  
  dependencies <- c("terra")
  
  install.packages(setdiff(dependencies, rownames(installed.packages())))
  
  # Validate parameters --------------------------------------------------------
  
  # Validate probability raster file
  if (!file.exists(probRasterFile))
    stop(paste0("Probability raster file not found: '", probRasterFile, "'."))
  
  # Create proportion raster ---------------------------------------------------
  
  # Load the probability raster
  probRaster <- terra::rast(probRasterFile)
  
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

  # Save proportion raster
  terra::writeRaster(propRaster, outputFile)
  
}

# Entrypoint for ArcGIS script tool
tool_exec <- function(in_params, out_params) {
  
  # if connected to internet, install latest version from github. 
  # This step will be skipped if commit hash matches
  
  # Install from current directory if version doesn't match.
  # ArcGIS runs R session in directory of current file
  if (!"LandslideSusceptibilityTW" %in% rownames(installed.packages()) |
      desc::desc_get_version("../DESCRIPTION") != 
      packageVersion("LandslideSusceptibilityTW")) {
    
    # Install dependencies
    
    # Install from source
    install.packages("..", 
                     repos = NULL, 
                     type = "source", 
                     dependencies = c("Imports", 
                                      "Suggests"))
  }
  
  createProportionRaster(
    probRasterFile = in_params[[1]],
    outputFile     = out_params[[1]] 
  )
  
  return(out_params)
}
