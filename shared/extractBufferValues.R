#' @title Extract buffer values
#' 
#' @description Extracts variable raster values within initiation and 
#' non-initiation buffers.
#' 
#' @param raster           A SpatRaster of explanatory variables
#' @param initBuffers      A SpatVector of initiation buffers
#' @param noninitBuffers   A SpatVector of non-initiation buffers
#' @param extractionMethod Method to use for extracting values from each buffer:
#'                         "all cells", "center cell", "max gradient cell", or 
#'                         "max plan cell"
#' 
#' @return A dataframe of extracted raster values with an additional "class" 
#' column.

extractBufferValues <- function(raster, initBuffers, noninitBuffers, extractionMethod) {
  
  # By default, extract all values from initiation and non-initiation buffers
  initValues <- terra::extract(raster, initBuffers)
  noninitValues <- terra::extract(raster, noninitBuffers)
  
  # Subset if a different buffer extraction method was requested
  if (extractionMethod == "center cell") {
    
    # Find center points of buffers
    initCenters <- terra::centroids(initBuffers)
    noninitCenters <- terra::centroids(noninitBuffers)
    
    # Get center point coordinates
    initCoords <- terra::geom(initCenters)[,c("x","y")]
    noninitCoords <- terra::geom(noninitCenters)[,c("x","y")]
    
    # Extract values from cells containing center points
    initValues <- terra::extract(raster, initCoords)
    noninitValues <- terra::extract(raster, noninitCoords)
    
  } else if (extractionMethod == "max gradient cell") {
    
    initValues <- aggregateBufferValues(initValues, "grad", max)
    noninitValues <- aggregateBufferValues(noninitValues, "grad", max)
    
  } else if (extractionMethod == "max plan cell") {
    
    initValues <- aggregateBufferValues(initValues, "plan", max)
    noninitValues <- aggregateBufferValues(noninitValues, "plan", max)
    
  }
  
  # Assign a classification value to each entry
  initValues$class <- rep("initiation", nrow(initValues))
  noninitValues$class <- rep("non-initiation", nrow(noninitValues))
  
  # Combine initiation and non-initiation entries into a single dataset
  dataset <- rbind(initValues, noninitValues)
  
  # Remove the "ID" column
  dataset$ID <- NULL
  
  # TODO: Remove cells that fall within initiation AND non-initiation buffers
  
  # Factor the classification variable values
  dataset$class <- factor(dataset$class)
  
  # Filter out entries with NA values
  dataset <- na.omit(dataset)
  
  return(dataset)
  
}


#' @title Aggregate buffer values
#' 
#' @description Groups values by buffer and--for each buffer--keeps the entry 
#' with the "fun" (min/max) variable value.
#' 
#' @param values  Values extracted from buffers, including an "ID" column
#' @param varName Name pattern of the variable to aggregate by
#' @param fun     A function to aggregate buffer values by: max, min
#' 
#' @return A dataframe of buffer values.

aggregateBufferValues <- function(values, varName, fun) {
  
  # Determine which variable to aggregate by based on provided name pattern
  valuesVarName <- names(values)[grepl(varName, names(values))][1]
  
  # Formula to group entries by buffer and return requested variable value
  fm <- as.formula(paste(valuesVarName, "~", "ID"))
  
  # For each buffer, keep the entry with the "fun" variable value
  aggregatedValues <- merge(
    aggregate(fm, max, data = values),
    values
  )
  
  return(aggregatedValues)
  
}