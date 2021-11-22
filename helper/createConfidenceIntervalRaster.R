#' @title Create a confidence interval raster
#' 
#' @description Creates a raster that stores the confidence interval bounds at 
#' each cell of a set of probability rasters.
#' 
#' @param probRaster A multi-layered raster of probability values from 0.0 to 1.0
#' @param p          Confidence interval p-value 

createConfidenceIntervalRaster <- function(x, p = 0.05) {
  n <- terra::nlyr(x)
  degreesFreedom <- n - 1
  tScore <- qt(p = p, df = degreesFreedom, lower.tail = FALSE)
  
  sd <- terra::app(x, fun = "sd") # Standard deviation
  se <- sd / sqrt(n)              # Standard error
  me <- se * tScore               # Margin of error
  
  mean <- terra::mean(x)          # Mean
  lower <- mean - me              # Interval lower bound
  upper <- mean + me              # Interval upper bound
  
  interval <- c(lower, upper)
  names(interval) <- c("lower", "upper")
  
  return(interval)
}