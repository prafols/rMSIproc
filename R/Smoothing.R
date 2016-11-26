#Smoothin methods


#' Smoothing.
#'
#' @param x the intensities of a a mass spectrum to smooth.
#' @param method the method to use for smoothing, available option are: "SavitzkyGolay".
#' @param ... specific parameters to each smoothing method
#'
#' Smooths a vector of data using the specified smoothing method.  
#'  
#' @return the smoothed data.
#' @export
#'
Smoothing <- function(x, method = "SavitzkyGolay", ...)
{
  y <- NULL
  if(method == "SavitzkyGolay")
  {
    y <- Smoothing_SavitzkyGolay(x, ...)
  }
  if(is.null(y))
  {
    stop("The selected method is not valid. Valid methods are: SavitzkyGolay")
  }
  return(y)
}