#' DetectPeaks'
#'
#' @param mass a vector containing the mass axis.
#' @param intensity a vector containinf the spectrum intensities.
#' @param SNR the minimum signal to noise ratio of retained peaks
#' @param WinSize the used windows size for peak detection
#'
#' @return a list containing mass, intensity and SNR fields of detected peaks.
#' @export
#'
DetectPeaks <- function(mass, intensity, SNR = 5, WinSize = 20)
{
  pm <- DetectPeaks_C(mass, intensity,SNR, WinSize)
  peaks <- list()
  peaks$mass <- pm["mass", ]
  peaks$intensity <- pm["intensity", ]
  peaks$SNR <- pm["SNR", ]
  return(peaks)
}