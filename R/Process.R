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
#' @param EnableSmoothing a boolean declaring if smoothing alignment must be performed.
#' @param SmoothingKernelSize size of smoothing kernel. NULL value may be specified if EnableSmoothing is FALSE.
#' @param EnableAlignment a boolean declaring if label-free alignment must be performed.
#' @param AlignmentIterations number of iterations of label-free alignment. The rMSI ramdisk will be overwritted with aligned data. NULL value may be specified if EnableAlignment is FALSE.
#' @param AlignmentMaxShiftppm the maximum shift that alignment can apply in ppm. NULL value may be specified if EnableAlignment is FALSE.
#' @param EnableCalibration a boolean declaring if mass calibration must be performed.
#' @param CalibrationPeakWin the windows size used for peak detection in calibration window.
#' @param EnablePeakPicking a boolean declaring if peak-pickin (and binning) must be performed.
#' @param SNR minimal singal to noise ratio of peaks to retain. NULL value may be specified if EnablePeakPicking is FALSE.
#' @param peakWindow windows size used for peak detection. Generally should be similar to peak with number of data points.  NULL value may be specified if EnablePeakPicking is FALSE.
#' @param peakUpSampling upsampling factor used in peak interpolation fo exact mass prediction. NULL value may be specified if EnablePeakPicking is FALSE.
#' @param UseBinning if true binned matrices are returned instead of peak lists.
#' @param BinTolerance the tolerance used to merge peaks to the same bin. It is recomanded to use the peak width in Da units. NULL value may be specified if EnablePeakPicking is FALSE.
#' @param BinFilter the peaks bins non detected in at least the BinFitler*TotalNumberOfPixels spectra will be deleted. NULL value may be specified if EnablePeakPicking is FALSE.
#' @param EnableSpectraNormalization if normalization must be applied.
#' @param EnableTICNorm if TIC normalization must be performed on spectra.
#' @param EnableMAXNorm if MAX normalization must be performed on spectra.
#' @param EnableTICAcqNorm if TICAcq normalization must be performed on spectra.
#' @param NumOfThreads the number number of threads used to process the data.
#' @param CalSpan the used span for the loess fittin used in mass calibration.
#' 
#'
#' @return  a named list containing:
#'             - The process image reference (procImg).
#'             - The results peak-picking (peakMat). This can be returned in two forms:
#'                 From1 (if binning is used) - a list containing three matrices (intensity, SNR and area) and a vector with a common mass axis.
#'                 Form2 (if NO binning is applied) - a list of detected peaks for each pixel.
#'             - The applied mass shifts in first alignment iteration (alignShifts).
#'
ProcessImage <- function(img, 
                         EnableSmoothing = T, SmoothingKernelSize = 5,
                         EnableAlignment = T, AlignmentIterations = 3, AlignmentMaxShiftppm = 200,
                         EnableCalibration = T, CalibrationPeakWin = 20,
                         EnablePeakPicking = T, SNR = 5, peakWindow = 10, peakUpSampling = 10, 
                         UseBinning = T, BinTolerance = 0.05, BinFilter = 0.05,
                         EnableSpectraNormalization = T, EnableTICNorm = T, EnableMAXNorm = T, EnableTICAcqNorm = T,
                         NumOfThreads = parallel::detectCores(), CalSpan = 0.75)
{
  pt <- Sys.time()
  
  dataInf <- getrMSIdataInfo(img)

  #Avoid using MALDIquant from here
  if(class(img$mean) == "MassSpectrum")
  {
    img$mean <- img$mean@intensity
  }
  
  #Apply Savitzky-Golay smoothing to RAW data and average spectrum
  if( EnableSmoothing )
  {
    cat("Running Savitzky-Golay Smoothing...\n")
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
  if( EnableAlignment )
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
  if( EnableSmoothing || EnableAlignment )
  {
    img$mean <- rMSI::AverageSpectrum(img)
  }
  
  #Do not count time while the user is in manual calibration window
  elap_1st_stage <- Sys.time() - pt
  
  #Manual calibration (user will be promp with calibration dialog)
  if( EnableCalibration )
  {
    img$mass <- CalibrationWindow( img$mass, img$mean, CalibrationPeakWin ,img$name, CalibrationSpan = CalSpan )
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
  if( EnablePeakPicking )
  {
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
    
    if(UseBinning)
    {
      cat("Replacing zero values in the binned peak matrix...\n")
      pkMatrix <- ReplacePeakMatrixZerosC(pkMatrix, 
                                          basePath = dataInf$basepath, 
                                          fileNames = dataInf$filenames, 
                                          mass = img$mass, 
                                          numRows = dataInf$nrows, 
                                          dataType = dataInf$datatype, 
                                          numOfThreads = NumOfThreads, 
                                          WinSize = peakWindow,  
                                          InterpolationUpSampling = peakUpSampling )
    }
    
  }
  else
  {
    pkMatrix <- NULL
  }
  
  #Calculate some normalizations
  if( EnableSpectraNormalization )
  {
    if( EnableTICNorm )
    {
      img <- rMSI::NormalizeTIC(img, remove_empty_pixels = T)
    }
    if( EnableMAXNorm )
    {
      img <- rMSI::NormalizeMAX(img, remove_empty_pixels = T)
    }
    if( EnableTICAcqNorm )
    {
      img <- rMSI::NormalizeByAcqDegradation(img)
    }
  }
  
  #Add a copy of img$pos to pkMatrix
  if( EnablePeakPicking && UseBinning)
  {
    pkMatrix <- FormatPeakMatrix(pkMatrix, img$pos,  c(nrow(img$pos)), c(img$name))
  }
  
  elap_2nd_stage <- Sys.time() - pt
  cat("Total used processing time:\n")
  print(elap_1st_stage + elap_2nd_stage)
  
  return ( list( procImg = img,   peakMat = pkMatrix, alignShifts = alngLags ))
}

#' FormatPeakMatrix.
#' Formats a C style peak matrix generated by MTPeakPicking::BinPeaks() to a rMSIprocPeakMatrix.
#'
#' @param cPeakMatrix a peak matrix with the same format as retured by MTPeakPicking::BinPeaks().
#' @param posMat a rMSI image pos matrix.
#' @param numPixels a vector including the number of pixels of each sample.
#' @param names a vector of strings with the name of each sample.
#'
#' @return the formated matrix.
#'
FormatPeakMatrix <- function (cPeakMatrix, posMat, numPixels, names)
{
  cPeakMatrix$pos <- posMat
  cPeakMatrix$numPixels <- numPixels
  cPeakMatrix$names <- names
  class(cPeakMatrix) <- "rMSIprocPeakMatrix"
  return(cPeakMatrix)
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
#'
MergePeakMatrices <- function( PeakMatrixList, binningTolerance = 0.05, binningFilter = 0.01, OffsetPosByX = F  )
{
  pt <- Sys.time()
  
  #Merge peak matrices
  pkMatrix <- MergePeakMatricesC( PeakMatrixList, binningTolerance, binningFilter )
  
  #Testing if normalizations arrays are concatenable (otherwise raise an error)
  normNames1 <- NULL
  if( !is.null(PeakMatrixList[[1]]$normalizations))
  {
    normNames1 <- names(PeakMatrixList[[1]]$normalizations)
  }
  for( i in 1:length(PeakMatrixList))
  {
    if( is.null(normNames1) != is.null(PeakMatrixList[[i]]$normalizations))
    {
      stop("Error: Normalizations is not present for all peak matrices. Please, provide peak matrices with same normalizations.\n")
    }
    if( !is.null(PeakMatrixList[[i]]$normalizations))
    {
      normNamesI <- names(PeakMatrixList[[i]]$normalizations)
      if( length(normNamesI) != length(normNames1) )
      {
        stop("Error: Peak matrices contains different numbers of normalizations.\n")
      }
      for( j  in 1:length(normNames1))
      {
        if( normNamesI[j] != normNames1[j] )
        {
          stop("Error: Peak matrices contains normalizations with differents names.\n")
        }
      }
    }
  }
  
  #Concatenate pos arrays and normalizations arrays 
  numOfPixels <- sum(unlist(lapply(PeakMatrixList, function(x){ nrow(x$pos) })))
  pkMatrix$pos <- matrix(nrow = numOfPixels, ncol = 2 )
  if( !is.null(normNames1) )
  {
    pkMatrix$normalizations <- data.frame( matrix( NA, nrow = numOfPixels, ncol = length(normNames1)) )
    names(pkMatrix$normalizations) <- normNames1
  }
  
  colnames(pkMatrix$pos) <- c("x", "y")
  iStart <- 1 #Index of matrix row
  MaxXant <- 0
  MaxYant <- 0
  pkMatrix$numPixels <- c()
  sampleNames <- c()
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
    
    #Set normalizations
    if( !is.null(pkMatrix$normalizations) )
    {
      pkMatrix$normalizations[ (iStart:iStop) , ] <- PeakMatrixList[[i]]$normalizations[,]
    }
    iStart <- iStop + 1
    
    #Set matrix names
    if( is.null(PeakMatrixList[[i]]$names) )
    {
      cat(paste("Warning: No sample name found in matrix", i, "\n"))
      sampleNames <- c(sampleNames, paste("unamed_sample_", i, sep = ""))
    }
    else
    {
      sampleNames <- c(sampleNames, PeakMatrixList[[i]]$names )
    }
  }
  
  elap <- Sys.time() - pt
  cat("Total used processing time:\n")
  print(elap)
  
  return( FormatPeakMatrix(pkMatrix, pkMatrix$pos, pkMatrix$numPixels, sampleNames))
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
#' @param overwriteRamdisk if the current ramdisk must be overwrited.
#' @param calibrationSpan the used span in the loess fitting for mass calibration.
#'
ProcessWizard <- function( deleteRamdisk = T, overwriteRamdisk = F, calibrationSpan = 0.75 )
{
  #Get processing params using a GUI
  procParams <- ImportWizardGui()
  if(is.null(procParams))
  {
    cat("Processing aborted\n")
    return()
  }
  
  #Get number of images
  brukerXmlCounters <- rep(0, length(procParams$data$source$xmlpath))
  if( procParams$data$source$type == "xmass" )
  {
    for( i in 1:length(procParams$data$source$xmlpath))
    {
      brukerXmlCounters[i] <- rMSI:::CountImagesInBrukerXml(procParams$data$source$xmlpath[i])
    }
    NumOfImages <- sum(brukerXmlCounters)
  }
  else
  {
    NumOfImages <- length(procParams$data$source$datapath)
  }
  
  brukerCounter_ant <- 0
  selBrukerXML <- 1
  selBrukerClass <- 0
  for( i in 1:NumOfImages)
  {
    cat(paste("Working on image", i, "of", NumOfImages, "\n"))
    
    #Load each image
    if( procParams$data$source$type == "xmass" )
    {
      if (brukerXmlCounters[selBrukerXML] >= i - brukerCounter_ant)
      {
        selBrukerClass <- selBrukerClass + 1
      }
      else
      {
        brukerCounter_ant <- brukerXmlCounters[selBrukerXML] + brukerCounter_ant
        selBrukerXML <- selBrukerXML + 1
        selBrukerClass <- 1
      }
      mImg <- rMSI::importBrukerXmassImg( procParams$data$source$datapath, procParams$data$pixelsize, procParams$data$source$xmlpath[selBrukerXML], procParams$data$source$spectrumpath, selected_img = selBrukerClass )
    }
    else
    {
      mImg <- rMSI::LoadMsiData( procParams$data$source$datapath[i], ff_overwrite = overwriteRamdisk )  
    }
    
    #Process data
    procData <- ProcessImage(img = mImg,
                             EnableSmoothing = procParams$smoothing$enabled, SmoothingKernelSize = procParams$smoothing$sgkernsize,
                             EnableAlignment = procParams$alignment$enabled, AlignmentIterations = procParams$alignment$iterations, AlignmentMaxShiftppm = procParams$alignment$maxshift,
                             EnableCalibration = procParams$calibration$enabled, CalibrationPeakWin = procParams$calibration$winsize,
                             EnablePeakPicking = procParams$peakpicking$enabled, SNR = procParams$peakpicking$snr, peakWindow = procParams$peakpicking$winsize, peakUpSampling = procParams$peakpicking$oversample,
                             UseBinning = T, BinTolerance = procParams$peakpicking$bintolerance, BinFilter = procParams$peakpicking$binfilter,
                             EnableSpectraNormalization = procParams$spectraNormalization$enabled, EnableTICNorm = procParams$spectraNormalization$TIC, EnableMAXNorm = procParams$spectraNormalization$MAX, EnableTICAcqNorm = procParams$spectraNormalization$AcqTIC,
                             NumOfThreads = procParams$nthreads, CalSpan = calibrationSpan )

    #Store MS image to a tar file
    if( procParams$smoothing$enabled || procParams$alignment$enabled || procParams$calibration$enabled || procParams$spectraNormalization$enabled || procParams$data$source$type != "tar"  )
    {
      rm(mImg) #Now all is stored in procData
      gc()
      imgName <- sub('\\..*$', '',procData$procImg$name) #Remove extensions of image name
      rMSI::SaveMsiData( procData$procImg, file.path(procParams$data$outpath,  paste(imgName,"-proc.tar", sep ="")))
    }
    else
    {
      imgName <- sub('\\..*$', '',mImg$name) #Remove extensions of image name
    }

    #Store peak matrix
    if( procParams$peakpicking$enabled )
    {
      #Append normalizations to
      if( !is.null( procData$procImg$normalizations) )
      {
        procData$peakMat$normalizations <-  as.data.frame(procData$procImg$normalizations )
      }

      StorePeakMatrix( file.path(procParams$data$outpath,  paste(imgName,"-peaks.zip", sep ="")), procData$peakMat)
    }

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
  
  writeLines(paste("Smoothing enabled = ", procParams$smoothing$enabled, sep ="" ), con = fObj)
  if(procParams$smoothing$enabled)
  {
    writeLines(paste("Smoothing SG kernel size = ", procParams$smoothing$sgkernsize, sep ="" ), con = fObj)
  }
  
  writeLines(paste("Alignment enabled = ", procParams$alignment$enabled, sep ="" ), con = fObj)
  if(procParams$alignment$enabled)
  {
    writeLines(paste("Alignment iterations = ", procParams$alignment$iterations, sep ="" ), con = fObj)
    writeLines(paste("Alignment max shift [ppm] = ", procParams$alignment$maxshift, sep ="" ), con = fObj)
  }
  
  writeLines(paste("Calibration enabled = ", procParams$calibration$enabled, sep ="" ), con = fObj)
  if(procParams$calibration$enabled)
  {
    #Reserverd for future expansion...
  }

  if(procParams$spectraNormalization$enabled)
  {
    writeLines(paste("Spectra TIC normalization enabled = ", procParams$spectraNormalization$TIC, sep ="" ), con = fObj)
    writeLines(paste("Spectra MAX normalization enabled = ", procParams$spectraNormalization$MAX, sep ="" ), con = fObj)
    writeLines(paste("Spectra TICAcq normalization enabled = ", procParams$spectraNormalization$AcqTIC, sep ="" ), con = fObj)
  }
  
  writeLines(paste("Peak-picking enabled = ", procParams$peakpicking$enabled, sep ="" ), con = fObj)
  if( procParams$peakpicking$enabled)
  {
    writeLines(paste("Peak-picking SNR threshold = ", procParams$peakpicking$snr, sep ="" ), con = fObj)
    writeLines(paste("Peak-picking detector window = ", procParams$peakpicking$winsize, sep ="" ), con = fObj)
    writeLines(paste("Peak-picking oversampling = ", procParams$peakpicking$oversample, sep ="" ), con = fObj)
        writeLines(paste("Peak-picking binning tolerance = ", procParams$peakpicking$bintolerance, sep ="" ), con = fObj)
    writeLines(paste("Peak-picking binning filter = ", procParams$peakpicking$binfilter, sep ="" ), con = fObj)
  }
  close(fObj)
}
