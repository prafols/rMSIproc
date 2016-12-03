#########################################################################
#     rMSIproc - R package for MSI data processing
#     Copyright (C) 2014 Pere Rafols Soler
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

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