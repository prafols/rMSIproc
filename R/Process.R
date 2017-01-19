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
#' @param AlignmentMaxShiftppm the maximum shift that alignment can apply in ppm.
#' @param SNR minimal singal to noise ratio of peaks to retain.
#' @param peakWindow windows size used for peak detection. Generally should be similar to peak with number of data points.
#' @param peakUpSampling upsampling factor used in peak interpolation fo exact mass prediction.
#' @param SmoothingKernelSize size of smoothing kernel, if zero smoothing is disabled.
#' @param UseManualCalibration if the manual calibration windows must be used.
#' @param UseBinning if true binned matrices are returned instead of peak lists.
#' @param BinTolerance the tolerance used to merge peaks to the same bin. It is recomanded to use the peak width in Da units.
#' @param BinFilter the peaks bins non detected in at least the BinFitler*TotalNumberOfPixels spectra will be deleted.
#' @param NumOfThreads the number number of threads used to process the data.
#' 
#'
#' @return  a named list containing:
#'             - The process image reference (procImg).
#'             - The results peak-picking (peakMat). This can be returned in two forms:
#'                 From1 (if binning is used) - a list containing three matrices (intensity, SNR and area) and a vector with a common mass axis.
#'                 Form2 (if NO binning is applied) - a list of detected peaks for each pixel.
#'             - The applied mass shifts in first alignment iteration (alignShifts).
#' 
#' @export
#'
ProcessImage <- function(img, AlignmentIterations = 0, AlignmentMaxShiftppm = 200, SNR = 5, peakWindow = 10, peakUpSampling = 10, 
                         SmoothingKernelSize = 5, UseManualCalibration = T,
                         UseBinning = T, BinTolerance = 0.05, BinFilter = 0.05, 
                         NumOfThreads = parallel::detectCores())
{
  pt <- Sys.time()
  
  dataInf <- getrMSIdataInfo(img)

  #Avoid using MALDIquant from here
  if(class(img$mean) == "MassSpectrum")
  {
    img$mean <- img$mean@intensity
  }
  
  #Apply Savitzky-Golay smoothing to RAW data and average spectrum
  if(SmoothingKernelSize > 0)
  {
    cat("Running Savitzkt-Golay Smoothing...\n")
    img$mean <- Smoothing_SavitzkyGolay(img$mean, SmoothingKernelSize)
    FullImageSmoothing(basePath = dataInf$basepath, 
                       fileNames = dataInf$filenames, 
                       massChannels = length(img$mass), 
                       numRows = dataInf$nrows,
                       dataType = dataInf$datatype, 
                       numOfThreads = NumOfThreads, 
                       SmoothingKernelSize = SmoothingKernelSize)
  }
  
  #Label-free Alignment
  if(AlignmentIterations > 0)
  {
    #Calculate reference spectrum for label free alignment
    refSpc <- InternalReferenceSpectrum(img)
    
    cat("Running Label-Free Alignment...\n")
    alngLags <- FullImageAlign(basePath = dataInf$basepath,
                               fileNames = dataInf$filenames, 
                               refSpectrum = refSpc, 
                               numRows = dataInf$nrows,
                               dataType = dataInf$datatype, 
                               numOfThreads = NumOfThreads, 
                               AlignmentIterations = AlignmentIterations,
                               AlignmentMaxShiftPpm = AlignmentMaxShiftppm)
  }
  else
  {
    alngLags <- NULL
  }
  
  #Recalculate mean spectrum
  img$mean <- rMSI::AverageSpectrum(img)
  
  #Do not count time while the user is in manual calibration window
  elap_1st_stage <- Sys.time() - pt
  
  #Manual calibration (user will be promp with calibration dialog)
  if( UseManualCalibration )
  {
    img$mass <- CalibrationWindow( img$mass, img$mean, img$name )
    if(is.null(img$mass))
    {
      rMSI::DeleteRamdisk(img)
      gc()
      stop("Aborted by user\n")
    }
  }
  
  #Reset elapset time counter
  pt <- Sys.time()
  
  #Peak-Picking and binning
  cat("Running Peak Picking...\n")
  pkMatrix <- FullImagePeakPicking(basePath = dataInf$basepath,
                                   fileNames = dataInf$filenames,
                                   mass = img$mass,  
                                   numRows = dataInf$nrows,
                                   dataType = dataInf$datatype, 
                                   numOfThreads = NumOfThreads, 
                                   SNR = SNR, 
                                   WinSize = peakWindow, 
                                   InterpolationUpSampling = peakUpSampling, 
                                   doBinning = UseBinning, 
                                   binningTolerance = BinTolerance, 
                                   binningFilter = BinFilter)
  
  #Calculate some normalizations
  img <- rMSI::NormalizeTIC(img, remove_empty_pixels = T)
  img <- rMSI::NormalizeMAX(img, remove_empty_pixels = T)
  img <- rMSI::NormalizeByAcqDegradation(img)
  
  #Add a copy of img$pos to pkMatrix
  if(UseBinning)
  {
    pkMatrix$pos <- img$pos
    pkMatrix$numPixels <- nrow(img$pos)
  }
  
  elap_2nd_stage <- Sys.time() - pt
  cat("Total used processing time:\n")
  print(elap_1st_stage + elap_2nd_stage)
  
  return ( list( procImg = img,   peakMat = pkMatrix, alignShifts = alngLags ))
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

#' ProcessWizard.
#' 
#' Imports and process MSI data using a friendly GUI.
#' Various images can be loaded and processed with a single execution.
#' Data can be in XMASS, tar (rMSI) or imzML format.
#' Processed data will be saved in a user specified directory.
#' The applied processing consists in:
#'   - Label-free aligment (various iterations can be performed, zero iterations means no alignment).
#'   - Peak-picking.
#'   - Peak-binning.
#'   - Mass calibration with internal reference compounds.
#' Processed data includes:
#'   - a .tar file with the processed data.
#'   - a rMSIproc formated matrices with binned peaks.
#'   - a plain text file with used processing parameters.
#' @param deleteRamdisk if the used ramdisks for MS images must be deleted for each image processing (will be deleted after saving it to .tar file).
#'  
#' @export
#'
ProcessWizard <- function( deleteRamdisk = T )
{
  #Get processing params using a GUI
  procParams <- ImportWizardGui()
  if(is.null(procParams))
  {
    cat("Processing aborted\n")
    return()
  }
  
  #Get number of images
  if( procParams$data$source$type == "xmass" )
  {
    NumOfImages <- length(procParams$data$source$xmlpath)
  }
  else
  {
    NumOfImages <- length(procParams$data$source$datapath)
  }
  
  for( i in 1:NumOfImages)
  {
    cat(paste("Working on image", i, "of", NumOfImages, "\n"))
    
    #Load each image
    if( procParams$data$source$type == "xmass" )
    {
      mImg <- rMSI::importBrukerXmassImg( procParams$data$source$datapath, procParams$data$pixelsize, procParams$data$source$xmlpath[i], procParams$data$source$spectrumpath )
    }
    else
    {
      mImg <- rMSI::LoadMsiData( procParams$data$source$datapath[i], ff_overwrite = T )  
    }
    
    #Process data
    procData <- ProcessImage(img = mImg, 
                 AlignmentIterations = procParams$alignment$iterations, AlignmentMaxShiftppm = procParams$alignment$maxshift,
                 SNR = procParams$peakpicking$snr, peakWindow = procParams$peakpicking$winsize, peakUpSampling = procParams$peakpicking$oversample, SmoothingKernelSize = procParams$smoothing$sgkernsize, 
                 UseBinning = T, BinTolerance = procParams$binning$tolerance, BinFilter = procParams$binning$filter, NumOfThreads = procParams$nthreads )
    
    rm(mImg) #Now all is stored in procData
    gc()
    
    #Store MS image to a tar file
    imgName <- sub('\\..*$', '',procData$procImg$name) #Remove extensions of image name
    rMSI::SaveMsiData( procData$procImg, file.path(procParams$data$outpath,  paste(imgName,"-proc.tar", sep ="")))
    
    #Store peak matrix
    StorePeakMatrix( file.path(procParams$data$outpath,  paste(imgName,"-peaks.zip", sep ="")), procData$peakMat)
    
    #Delete data and clear memory
    if(deleteRamdisk)
    {
      rMSI::DeleteRamdisk(procData$procImg)
    }
    gc()
  }
  
  #Store processing params
  SaveProcessingParams( procParams, file.path(procParams$data$outpath, "processing_parameters.txt"  ))
  
  cat("All images were successfully processed\n")
}

#' SaveProcessingParams.
#' 
#' Save all parameters in a list of processing params generated using ImportWizardGui() function.
#' Parameters will be saved in a plain text file.
#'
#' @param procParams a list of parameters.
#' @param filepath a full path where params will be stored
#'
#' @export
#'
SaveProcessingParams <- function( procParams, filepath)
{
  fObj <- file(description = filepath, open = "w" )
  writeLines("Processing parameters", con = fObj)
  writeLines(paste("Data source type = ", procParams$data$source$type, sep ="" ), con = fObj)
  if( procParams$data$source$type == "xmass" )
  {
    writeLines(paste("Data source directory = ", procParams$data$source$datapath, sep ="" ), con = fObj)
    for(i in 1:length( procParams$data$source$xmlpath))
    {
      writeLines(paste("Data source XML file", i," = ", procParams$data$source$xmlpath[i], sep ="" ), con = fObj)
    }
    writeLines(paste("Data source ref spectrum = ", procParams$data$source$spectrumpath, sep ="" ), con = fObj)
  }
  else
  {
    for(i in 1:length( procParams$data$source$datapath))
    {
      writeLines(paste("Data source file", i," = ", procParams$data$source$datapath[i], sep ="" ), con = fObj)
    }
  }
  writeLines(paste("Data output directory = ", procParams$data$outpath, sep ="" ), con = fObj)
  writeLines(paste("Alignment iterations = ", procParams$alignment$iterations, sep ="" ), con = fObj)
  writeLines(paste("Alignment max shift [ppm] = ", procParams$alignment$maxshift, sep ="" ), con = fObj)
  
  writeLines(paste("Peak-picking SNR threshold = ", procParams$peakpicking$snr, sep ="" ), con = fObj)
  writeLines(paste("Peak-picking detector window = ", procParams$peakpicking$winsize, sep ="" ), con = fObj)
  writeLines(paste("Peak-picking oversampling = ", procParams$peakpicking$oversample, sep ="" ), con = fObj)
  writeLines(paste("Peak-picking SG kernel = ", procParams$peakpicking$sgkernsize, sep ="" ), con = fObj)
  
  writeLines(paste("Peak-binning tolerance = ", procParams$binning$tolerance, sep ="" ), con = fObj)
  writeLines(paste("Peak-binning filter = ", procParams$binning$filter, sep ="" ), con = fObj)
  close(fObj)
}
