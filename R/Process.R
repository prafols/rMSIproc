
#' ProcessImage.
#' 
#' Perform all image pre-processing using a multi-threading implementation.
#' If aligment is used then the hdd files are overwirted with aligned data.
#' 
#' @param img an rMSI data object to process.
#' @param performAlignment if is true FFT based aligment will be used and the ramdisk overwritted with alignment data.
#' @param SNR minimal singal to noise ratio of peaks to retain.
#' @param peakWindow windows size used for peak detection. Generally should be similar to peak with number of data points.
#' @param peakUpSampling upsampling factor used in peak interpolation fo exact mass prediction.
#' @param SmoothingKernelSize size of smoothing kernel.
#' @param BinTolerance the tolerance used to merge peaks to the same bin. It is recomanded to use the peak width in Da units.
#' @param BinFilter the peaks bins non detected in at least the BinFitler*TotalNumberOfPixels spectra will be deleted.
#' @param NumOfThreads the number number of threads used to process the data.
#' 
#'
#' @return a intensity matrix where each row corresponds to an spectrum.
#' @export
#'
ProcessImage <- function(img, performAlignment = F, SNR = 5, peakWindow = 10, peakUpSampling = 10, 
                         SmoothingKernelSize = 5, BinTolerance = 0.05, BinFilter = 0.1, NumOfThreads = parallel::detectCores())
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
                               NumOfThreads, performAlignment, SNR, peakWindow, peakUpSampling, SmoothingKernelSize, BinTolerance, BinFilter)
  
  elap <- Sys.time() - pt
  cat("Total used processing time:\n")
  print(elap)
  
  return (pkMatrix )
  
}