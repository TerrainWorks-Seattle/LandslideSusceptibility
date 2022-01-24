#' @title Summarize a data frame
#' 
#' @details Calculates summary statistics for each data frame column and stores
#' them in a formatted matrix.
#' 
#' @param x the \code{data.frame} to summarize
#' 
#' @return a \code{matrix} of summary statistics
#' @export

summarizeDataFrame <- function(x) {
  
  # Calculate summary statistics for each column
  min <- as.data.frame(lapply(x, min, na.rm = TRUE))
  max <- as.data.frame(lapply(x, max, na.rm = TRUE))
  range <- as.data.frame(lapply(x, function(y) {
    max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
  }))
  mean <- as.data.frame(lapply(x, mean, na.rm = TRUE))
  sd <- as.data.frame(lapply(x, sd, na.rm = TRUE))
  
  # Create summary table
  df <- rbind(min, max, range, mean, sd)
  df <- t(df)
  summary <- as.matrix(df)
  colnames(summary) <- c("min", "max", "range", "mean", "sd")
  
  return(summary)
  
}