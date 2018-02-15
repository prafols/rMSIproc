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

#' AverageSpectrum.
#' 
#' Calculates the dataset average spectrum.
#' The average spectrum is the spectrum produced by weighting all the dataset spectra.
#'
#' @param img An rMSI objects.
#' @param NumOfThreads Number of threads.The Default value is the number of cores of the machine. 
#'
#' @return The average spectrum.
#' 
MTAverageSpectrum <- function(img, NumOfThreads = parallel::detectCores())
{
  dataInf <- getrMSIdataInfo(img) #Info from the img

  AvgSptr <- NULL
  AvgSptr<-AverageSpectrumC(dataInf$filenames, dataInf$masschannels, dataInf$nrows, dataInf$datatype, NumOfThreads)
  
  return(AvgSptr)
}
