# Create a dataset containing variable values and classes of cells within 
# initiation and non-initiation buffers.

extractBufferValues <- function(
  varsRaster,            # A SpatRaster of explanatory variables
  initiationBuffers,     # A SpatVector of initiation buffers
  noninitiationBuffers,  # A SpatVector of non-initiation buffers,
  bufferExtractionMethod # Method to use for extracting buffer values 
) {
  # By default, extract all values from initiation and non-initiation buffers
  initiationValues <- terra::extract(varsRaster, initiationBuffers)
  noninitiationValues <- terra::extract(varsRaster, noninitiationBuffers)
  
  # Subset if a different buffer extraction method was requested
  if (bufferExtractionMethod == "center cell") {
    # Find center points of buffers
    initiationCenters <- terra::centroids(initiationBuffers)
    noninitiationCenters <- terra::centroids(noninitiationBuffers)
    
    # Get center point coordinates
    initiationCoords <- terra::geom(initiationCenters)[,c("x","y")]
    noninitiationCoords <- terra::geom(noninitiationCenters)[,c("x","y")]
    
    # Extract values from cells containing center points
    initiationValues <- terra::extract(varsRaster, initiationCoords)
    noninitiationValues <- terra::extract(varsRaster, noninitiationCoords)
  } else if (bufferExtractionMethod == "max gradient cell") {
    initiationValues <- aggregateBufferValues(initiationValues, "grad", max)
    noninitiationValues <- aggregateBufferValues(noninitiationValues, "grad", max)
  } else if (bufferExtractionMethod == "max plan cell") {
    initiationValues <- aggregateBufferValues(initiationValues, "plan", max)
    noninitiationValues <- aggregateBufferValues(noninitiationValues, "plan", max)
  }
  
  # Assign a classification value to each entry
  initiationValues$class <- rep("initiation", nrow(initiationValues))
  noninitiationValues$class <- rep("non-initiation", nrow(noninitiationValues))
  
  # Combine initiation and non-initiation entries into a single dataset
  dataset <- rbind(initiationValues, noninitiationValues)
  
  # Remove the "ID" column
  dataset$ID <- NULL
  
  # Factor the classification variable values
  dataset$class <- factor(dataset$class)
  
  # Filter out entries with NA values
  dataset <- na.omit(dataset)
  
  return(dataset)
}

# Groups values by buffer and--for each buffer--keeps the entry with the "fun"
# (min/max) variable value.

aggregateBufferValues <- function(
  bufferValues, # Values extracted from buffers, including an "ID" column
  varName,      # Name pattern of the variable to aggregate by
  fun           # Function to aggregate buffer values by
) {
  # Determine which variable to aggregate by based on provided name pattern
  gradientVarName <- names(bufferValues)[grepl("grad", names(bufferValues))][1]
  
  # Formula to group entries by buffer and return requested variable value
  fm <- as.formula(paste(gradientVarName, "~", "ID"))
  
  # For each buffer, keep the entry with the "fun" variable value
  aggregatedValues <- merge(
    aggregate(fm, max, data = bufferValues),
    bufferValues
  )
  
  return(aggregatedValues)
}