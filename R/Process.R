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

#' ProcessImage.
#' 
#' Perform all image pre-processing using a multi-threading implementation.
#' If aligment is used then the hdd files are overwirted with aligned data.
#' A recomeneded value of aligment ietarations is 3.
#' 
#' @param img an rMSI data object to process.
#' @param AlignmentIterations if grather than zero FFT based aligment will be used and the ramdisk overwritted with alignment data.
#' @param SNR minimal singal to noise ratio of peaks to retain.
#' @param peakWindow windows size used for peak detection. Generally should be similar to peak with number of data points.
#' @param peakUpSampling upsampling factor used in peak interpolation fo exact mass prediction.
#' @param SmoothingKernelSize size of smoothing kernel.
#' @param UseBinning if true binned matrices are returned instead of peak lists.
#' @param BinTolerance the tolerance used to merge peaks to the same bin. It is recomanded to use the peak width in Da units.
#' @param BinFilter the peaks bins non detected in at least the BinFitler*TotalNumberOfPixels spectra will be deleted.
#' @param NumOfThreads the number number of threads used to process the data.
#' 
#'
#' @return a intensity matrix where each row corresponds to an spectrum.
#' @export
#'
ProcessImage <- function(img, AlignmentIterations = 0, SNR = 5, peakWindow = 10, peakUpSampling = 10, 
                         SmoothingKernelSize = 5, 
                         UseBinning = T, BinTolerance = 0.05, BinFilter = 0.05, 
                         NumOfThreads = parallel::detectCores())
{
  pt <- Sys.time()
  
  dataInf <- getrMSIdataInfo(img)

  
  #TODO add normalization and calibration here!
  
  if(class(img$mean) == "MassSpectrum")
  {
    refSpc <- img$mean@intensity
  }
  else
  {
    refSpc <- img$mean
  }
  
  pkMatrix <- FullImageProcess(dataInf$basepath, dataInf$filenames, img$mass, refSpc, dataInf$nrows, dataInf$datatype,
                               NumOfThreads, AlignmentIterations, SNR, peakWindow, peakUpSampling,
                               SmoothingKernelSize, UseBinning, BinTolerance, BinFilter)
  
  #Add a copy of img$pos to pkMatrix
  if(UseBinning)
  {
    pkMatrix$pos <- img$pos
    pkMatrix$numPixels <- nrow(img$pos)
  }
  
  elap <- Sys.time() - pt
  cat("Total used processing time:\n")
  print(elap)
  
  return (pkMatrix )
}

#' MergePeakMatrices.
#' 
#' Merges a list containing various peak matrices in a single peak matrix.
#' The rMSIproc binning method is used to calculate the new masses.
#'
#' @param PeakMatrixList A list of various peak matrix objexts produced using rMSIproc.
#' @param binningTolerance the tolerance used to merge peaks to the same bin. It is recomanded to use the peak width in Da units.
#' @param binningFilter the peaks bins non detected in at least the BinFitler*TotalNumberOfPixels spectra will be deleted.
#' @param OffsetPosByX if true the pos matrices are concatenated by offseting in X direction, if false Y direction is used.
#'
#' @return a intensity matrix where each row corresponds to an spectrum.
#' @export
#'
MergePeakMatrices <- function( PeakMatrixList, binningTolerance = 0.05, binningFilter = 0.01, OffsetPosByX = F  )
{
  pt <- Sys.time()
  
  #Merge peak matrices
  pkMatrix <- MergePeakMatricesC( PeakMatrixList, binningTolerance, binningFilter )
  
  #Concatenate pos arrays
  numOfPixels <- sum(unlist(lapply(PeakMatrixList, function(x){ nrow(x$pos) })))
  pkMatrix$pos <- matrix(nrow = numOfPixels, ncol = 2 )
  colnames(pkMatrix$pos) <- c("x", "y")
  iStart <- 1 #Index of matrix row
  MaxXant <- 0
  MaxYant <- 0
  pkMatrix$numPixels <- c()
  for( i in 1:length(PeakMatrixList))
  {
    iStop <- iStart + nrow(PeakMatrixList[[i]]$pos) - 1
    pkMatrix$pos[ (iStart:iStop ), ] <- PeakMatrixList[[i]]$pos
    pkMatrix$numPixels <- c(pkMatrix$numPixels, nrow(PeakMatrixList[[i]]$pos))
    if(OffsetPosByX)
    {
      pkMatrix$pos[ (iStart:iStop ), "x"] <- pkMatrix$pos[ (iStart:iStop ), "x"] + MaxXant
    }
    else
    {
      pkMatrix$pos[ (iStart:iStop ), "y"] <- pkMatrix$pos[ (iStart:iStop ), "y"] + MaxYant
    }
    MaxXant <- max( pkMatrix$pos[ (iStart:iStop ), "x"]  )
    MaxYant <- max( pkMatrix$pos[ (iStart:iStop ), "y"]  )
    iStart <- iStop + 1
  }
  
  elap <- Sys.time() - pt
  cat("Total used processing time:\n")
  print(elap)
  
  return( pkMatrix)
}
