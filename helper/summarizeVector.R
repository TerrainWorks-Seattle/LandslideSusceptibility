#' @title Summarize a numeric vector
#' 
#' @details Calculates summary statistics for a numeric vector and stores
#' them in a formatted matrix.
#' 
#' @param x the numeric vector to summarize
#' 
#' @return a \code{matrix} of summary statistics

summarizeVector <- function(x) {
  
  if (!is.numeric(x)) return(NULL)
    
  # Calculate summary statistics
  min <- min(x, na.rm = TRUE)
  max <- max(x, na.rm = TRUE)
  range <- diff(c(min, max))
  mean <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  
  # Create summary table
  summary <- matrix(c(min, max, range, mean, sd), ncol = 1)
  summary <- t(summary)
  colnames(summary) <- c("min", "max", "range", "mean", "sd")
  
  return(summary)
  
}