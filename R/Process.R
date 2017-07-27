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
#' @param img an rMSI data object to process or a list of rMSI objects if various datasets must merged for processing.
#' @param EnableSmoothing a boolean declaring if smoothing alignment must be performed.
#' @param SmoothingKernelSize size of smoothing kernel. NULL value may be specified if EnableSmoothing is FALSE.
#' @param EnableAlignment a boolean declaring if label-free alignment must be performed.
#' @param AlignmentIterations number of iterations of label-free alignment. The rMSI ramdisk will be overwritted with aligned data. NULL value may be specified if EnableAlignment is FALSE.
#' @param AlignmentMaxShiftppm the maximum shift that alignment can apply in ppm. NULL value may be specified if EnableAlignment is FALSE.
#' @param AlignmentBilinear if TRUE the biliniar algiment mode will be used insetad of linear.
#' @param AlignmentRefLow the relative location of the spectrum where the bottom part correlation is calculated.
#' @param AlignmentRefMid the relative location of the spectrum where the central part correlation is calculated (only for bilinear mode).
#' @param AlignmentRefHigh the relative location of the spectrum where the top part correlation is calculated.
#' @param AlignmentOversampling the oversampling used in spectrum scale/shift to provide better accuracy.
#' @param EnableCalibration a boolean declaring if mass calibration must be performed.
#' @param CalibrationPeakWin the windows size used for peak detection in calibration window.
#' @param EnablePeakPicking a boolean declaring if peak-pickin (and binning) must be performed.
#' @param SNR minimal singal to noise ratio of peaks to retain. NULL value may be specified if EnablePeakPicking is FALSE.
#' @param peakWindow windows size used for peak detection. Generally should be similar to peak with number of data points.  NULL value may be specified if EnablePeakPicking is FALSE.
#' @param peakUpSampling upsampling factor used in peak interpolation fo exact mass prediction. NULL value may be specified if EnablePeakPicking is FALSE.
#' @param UseBinning if true binned matrices are returned instead of peak lists.
#' @param BinTolerance the tolerance used to merge peaks to the same bin. It is recomanded to use the half of peak width in ppm units. NULL value may be specified if EnablePeakPicking is FALSE.
#' @param BinFilter the peaks bins non detected in at least the BinFitler*TotalNumberOfPixels spectra will be deleted. NULL value may be specified if EnablePeakPicking is FALSE.
#' @param EnableSpectraNormalization if normalization must be applied.
#' @param EnableTICNorm if TIC normalization must be performed on spectra.
#' @param EnableRMSNorm if RMS normalization must be performed on spectra.
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
                         EnableAlignment = T, AlignmentIterations = 3, AlignmentMaxShiftppm = 200, AlignmentBilinear = F,
                         AlignmentRefLow = 0, AlignmentRefMid = 0.5, AlignmentRefHigh = 1, AlignmentOversampling = 2,
                         EnableCalibration = T, CalibrationPeakWin = 20,
                         EnablePeakPicking = T, SNR = 5, peakWindow = 10, peakUpSampling = 10, 
                         UseBinning = T, BinTolerance = 100, BinFilter = 0.05, BinToleranceUsingPPM = T,
                         EnableSpectraNormalization = T, EnableTICNorm = T, EnableRMSNorm = T, EnableMAXNorm = T, EnableTICAcqNorm = T,
                         NumOfThreads = parallel::detectCores(), CalSpan = 0.75)
{
  pt <- Sys.time()
  
  
  if(class(img) == "rMSIObj")
  {
    #Single image processing
    dataInf <- getrMSIdataInfo(img)
    img_list <- list(img)
    rm(img)
    bReturnAList <- FALSE #Ensure the return the same data structre as it was provided
  }
  else
  {
    #Multiple image processing
    img_list <- img
    rm(img)
    for(i in 1:length(img_list))
    {
      if( class(img_list[[i]]) != "rMSIObj" )
      {
        stop("Error: Found an object that is not class of rMSIObj, aboring...\n")
      }
    }
    dataInf <- getrMSIdataInfoMultipleDataSets(img_list)
    bReturnAList <- TRUE #Ensure the return the same data structre as it was provided
  }
  
  #Avoid using MALDIquant from here
  for(i in 1:length(img_list))
  {
    if(class(img_list[[i]]$mean) == "MassSpectrum")
    {
      img_list[[i]]$mean <- img_list[[i]]$mean@intensity
    }
  }
  
  #Apply Savitzky-Golay smoothing to RAW data and average spectrum
  if( EnableSmoothing )
  {
    cat("Running Savitzky-Golay Smoothing...\n")
    #The ff files must be closed befor running the Cpp code
    for( i in 1:length(img_list))
    {
      lapply(img_list[[i]]$data, function(x){ ff::close.ff(x) })
    }
    
    FullImageSmoothing(fileNames = dataInf$filenames, 
                       massChannels = dataInf$masschannels, 
                       numRows = dataInf$nrows,
                       dataType = dataInf$datatype, 
                       numOfThreads = NumOfThreads, 
                       SmoothingKernelSize = SmoothingKernelSize)
    
    #The ff file must be re-open to continue
    for( i in 1:length(img_list))
    {
      lapply(img_list[[i]]$data, function(x){ ff::open.ff(x) })
    }
  }
  
  #Label-free Alignment
  if( EnableAlignment )
  {
    #Calculate reference spectrum for label free alignment
    refSpc <- InternalReferenceSpectrumMultipleDatasets(img_list)
    
    cat("Running Label-Free Alignment...\n")
    #The ff file must be closed befor running the Cpp code
    for( i in 1:length(img_list))
    {
      lapply(img_list[[i]]$data, function(x){ ff::close.ff(x) })
    }
    
    alngLags <- FullImageAlign(fileNames = dataInf$filenames, 
                               refSpectrum = refSpc, 
                               numRows = dataInf$nrows,
                               dataType = dataInf$datatype, 
                               numOfThreads = NumOfThreads, 
                               AlignmentBilinear = AlignmentBilinear,
                               AlignmentIterations = AlignmentIterations,
                               AlignmentMaxShiftPpm = AlignmentMaxShiftppm,
                               RefLow = AlignmentRefLow,
                               RefMid = AlignmentRefMid, 
                               RefHigh = AlignmentRefHigh, 
                               OverSampling = AlignmentOversampling )
    
    #The ff file must be re-open to continue
    for( i in 1:length(img_list))
    {
      lapply(img_list[[i]]$data, function(x){ ff::open.ff(x) })
    }
  }
  else
  {
    alngLags <- NULL
  }
  
  #Recalculate mean spectrum
  if( EnableSmoothing || EnableAlignment )
  {
    for( i in 1:length(img_list))
    {
      img_list[[i]]$mean <- rMSI::AverageSpectrum(img_list[[i]]) 
    }
  }
  
  #Do not count time while the user is in manual calibration window
  elap_1st_stage <- Sys.time() - pt
  
  #Manual calibration (user will be promp with calibration dialog)
  if( EnableCalibration )
  {
    cal_intensity_spc <- rep(0, dataInf$masschannels) 
    for( i in 1:length(img_list))
    {
      cal_intensity_spc <- cal_intensity_spc + img_list[[i]]$mean
    }
    cal_intensity_spc <- cal_intensity_spc/length(img_list)
    
    if(length(img_list) == 1)
    {
      str_cal_title <- img_list[[1]]$name
    }
    else
    {
      str_cal_title <- paste("Merged data of", length(img_list), "datasets")
    }
    
    common_mass <- CalibrationWindow( img_list[[1]]$mass, cal_intensity_spc, CalibrationPeakWin , str_cal_title, CalibrationSpan = CalSpan )
    if(is.null(common_mass))
    {
      for( i in 1:length(img_list))
      {
        rMSI::DeleteRamdisk(img_list[[i]])
      }
      gc()
      stop("Aborted by user\n")
    }
    
    #Store the common mass axis for all images
    for( i in 1:length(img_list))
    {
      img_list[[i]]$mass <- common_mass
    }
  }
  
  #Reset elapset time counter
  pt <- Sys.time()
  
  #Peak-Picking and binning
  if( EnablePeakPicking )
  {
    cat("Running Peak Picking...\n")
    pkMatrix <- FullImagePeakPicking(fileNames = dataInf$filenames,
                                     mass =  img_list[[1]]$mass,  
                                     numRows = dataInf$nrows,
                                     dataType = dataInf$datatype, 
                                     numOfThreads = NumOfThreads, 
                                     SNR = SNR, 
                                     WinSize = peakWindow, 
                                     InterpolationUpSampling = peakUpSampling, 
                                     doBinning = UseBinning, 
                                     binningTolerance = BinTolerance, 
                                     binningFilter = BinFilter,
                                     binningIn_ppm = BinToleranceUsingPPM)
    
    if(UseBinning)
    {
      cat("Replacing zero values in the binned peak matrix...\n")
      pkMatrix <- ReplacePeakMatrixZeros(pkMatrix, 
                                          fileNames = dataInf$filenames, 
                                          mass = img_list[[1]]$mass, 
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
    mergedNormList <- list()
    if( EnableTICNorm )
    {
      mergedNormList$TIC <- c()
    }
    if( EnableRMSNorm )
    {
      mergedNormList$RMS <- c()
    }
    if( EnableMAXNorm )
    {
      mergedNormList$MAX <- c()
    }
    if( EnableTICAcqNorm )
    {
      mergedNormList$AcqTic <- c()
    }
    
    for( i in 1:length(img_list))
    {
      cat(paste0("Normalizations for image ", i, "/", length(img_list), "\n"))
      if( EnableTICNorm )
      {
        img_list[[i]] <- rMSI::NormalizeTIC(img_list[[i]], remove_empty_pixels = F)
        mergedNormList$TIC <- c(mergedNormList$TIC,  img_list[[i]]$normalizations$TIC)
      }
      if( EnableRMSNorm )
      {
        img_list[[i]] <- rMSI::NormalizeRMS(img_list[[i]], remove_empty_pixels = F)
        mergedNormList$RMS <- c(mergedNormList$RMS,  img_list[[i]]$normalizations$RMS)
      }
      if( EnableMAXNorm )
      {
        img_list[[i]] <- rMSI::NormalizeMAX(img_list[[i]], remove_empty_pixels = F)
        mergedNormList$MAX <- c(mergedNormList$MAX,  img_list[[i]]$normalizations$MAX)
      }
      if( EnableTICAcqNorm )
      {
        img_list[[i]] <- rMSI::NormalizeByAcqDegradation(img_list[[i]])
        mergedNormList$AcqTic <- c(mergedNormList$AcqTic,  img_list[[i]]$normalizations$AcqTic)
      }
    }
    
    #Append normalizations to the peak matrix
    if(!is.null(pkMatrix))
    {
      pkMatrix$normalizations <-  as.data.frame(mergedNormList)
    }
    rm(mergedNormList)
  }
  
  #Add a copy of img$pos to pkMatrix
  if( EnablePeakPicking && UseBinning)
  {
    mergedNames <- unlist(lapply(img_list, function(x){ return(x$name) }))
    mergedNumPixels <- unlist(lapply(img_list, function(x){ return(nrow(x$pos)) }))
    mergedPos <- matrix(ncol = 2, nrow = sum(mergedNumPixels))
    mergedMotors <- matrix(ncol = 2, nrow = sum(mergedNumPixels))
    mergedUUIDs <- unlist(lapply(img_list, function(x){ return(x$uuid) }))
    colnames(mergedPos) <- c("x", "y")
    colnames(mergedMotors) <- c("x", "y")
    istart <- 1
    for( i in 1:length(img_list))
    {
      istop <- istart + nrow(img_list[[i]]$pos) - 1
      mergedPos[ istart:istop , "x"] <- img_list[[i]]$pos[, "x"]
      mergedPos[ istart:istop , "y"] <- img_list[[i]]$pos[, "y"]
      
      if(!is.null(img_list[[i]]$posMotors))
      {
        mergedMotors[ istart:istop , "x"] <- img_list[[i]]$posMotors[, "x"]
        mergedMotors[ istart:istop , "y"] <- img_list[[i]]$posMotors[, "y"]
      }
      else
      {
        mergedMotors[ istart:istop , "x"] <- img_list[[i]]$pos[, "x"]
        mergedMotors[ istart:istop , "y"] <- img_list[[i]]$pos[, "y"]
      }
      istart <- istop + 1 
    }
    pkMatrix <- FormatPeakMatrix(pkMatrix, mergedPos,  mergedNumPixels, mergedNames, mergedUUIDs, mergedMotors) 
  }
  
  elap_2nd_stage <- Sys.time() - pt
  cat("Total used processing time:\n")
  print(elap_1st_stage + elap_2nd_stage)
  
  return_list <- list()
  if( bReturnAList )
  {
    return_list$procImg <- img_list
  }
  else
  {
    return_list$procImg <- img_list[[1]]
  }
  
  return_list$peakMat <- pkMatrix
  return_list$alignShifts <- alngLags
  return ( return_list )
}

#' FormatPeakMatrix.
#' Formats a C style peak matrix generated by MTPeakPicking::BinPeaks() to a rMSIprocPeakMatrix.
#' If the original motors matrix posMotors is not provided, a copy of posMat will be used.
#'
#' @param cPeakMatrix a peak matrix with the same format as retured by MTPeakPicking::BinPeaks().
#' @param posMat a rMSI image pos matrix.
#' @param numPixels a vector including the number of pixels of each sample.
#' @param names a vector of strings with the name of each sample.
#' @param uuid a vector of img UUID to be also stored in peak matrices
#' @param posMotors a rMSI image original motros coordinates matrix.
#'
#' @return the formated matrix.
#'
FormatPeakMatrix <- function (cPeakMatrix, posMat, numPixels, names, uuid, posMotors = NULL)
{
  cPeakMatrix$pos <- posMat
  cPeakMatrix$numPixels <- numPixels
  cPeakMatrix$names <- names
  cPeakMatrix$uuid <- uuid
  if(!is.null(posMotors))
  {
    cPeakMatrix$posMotors <- posMotors
  }
  else
  {
    cPeakMatrix$posMotors <- posMat
  }
  class(cPeakMatrix) <- "rMSIprocPeakMatrix"
  return(cPeakMatrix)
}

#' MergePeakMatrices.
#' 
#' Merges a list containing various peak matrices in a single peak matrix.
#' The rMSIproc binning method is used to calculate the new masses.
#'
#' @param PeakMatrixList A list of various peak matrix objexts produced using rMSIproc.
#' @param binningTolerance the tolerance used to merge peaks to the same bin. It is recomanded to use the half of the peak width in ppm units.
#' @param binningFilter the peaks bins non detected in at least the BinFitler*TotalNumberOfPixels spectra will be deleted.
#'
#' @return a intensity matrix where each row corresponds to an spectrum.
#'
MergePeakMatrices <- function( PeakMatrixList, binningTolerance = 100, binningFilter = 0.01  )
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
  mergedUUIDs <- unlist(lapply(PeakMatrixList, function(x){x$uuid}))
  numOfPixels <- sum(unlist(lapply(PeakMatrixList, function(x){ nrow(x$pos) })))
  pkMatrix$pos <- matrix(nrow = numOfPixels, ncol = 2 )
  pkMatrix$posMotors <- matrix(nrow = numOfPixels, ncol = 2 )
  if( !is.null(normNames1) )
  {
    pkMatrix$normalizations <- data.frame( matrix( NA, nrow = numOfPixels, ncol = length(normNames1)) )
    names(pkMatrix$normalizations) <- normNames1
  }
  
  colnames(pkMatrix$pos) <- c("x", "y")
  colnames(pkMatrix$posMotors) <- c("x", "y")
  iStart <- 1 #Index of matrix row
  pkMatrix$numPixels <- c()
  sampleNames <- c()
  for( i in 1:length(PeakMatrixList))
  {
    iStop <- iStart + nrow(PeakMatrixList[[i]]$pos) - 1
    pkMatrix$pos[ (iStart:iStop ), ] <- PeakMatrixList[[i]]$pos
    if(is.null(PeakMatrixList[[i]]$posMotors))
    {
      pkMatrix$posMotors[ (iStart:iStop ), ] <- PeakMatrixList[[i]]$pos
    }
    else
    {
      pkMatrix$posMotors[ (iStart:iStop ), ] <- PeakMatrixList[[i]]$posMotors
    }
    pkMatrix$numPixels <- c(pkMatrix$numPixels, nrow(PeakMatrixList[[i]]$pos))
    
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
  
  return( FormatPeakMatrix(pkMatrix, pkMatrix$pos, pkMatrix$numPixels, sampleNames, mergedUUIDs, pkMatrix$posMotors))
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
  
  #Depending on the dataset merge bit the whole processing will be executed differnt times
  if(procParams$mergedatasets)
  {
    NumOfProcRuns <- 1
    NumOfImages2Merge <- NumOfImages
  }
  else
  {
    NumOfProcRuns <- NumOfImages
    NumOfImages2Merge <- 1
  }
  
  #Create structured list with all of the data paths
  dataPaths <- list()
  brukerCounter_ant <- 0
  selBrukerXML <- 1
  selBrukerClass <- 0
  for( i in 1:NumOfImages)
  {
    dataPaths[[i]] <- list()
    
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
      dataPaths[[i]]$name <- paste(basename(procParams$data$source$xmlpath[selBrukerXML]), "_", rMSI::ParseBrukerXML(procParams$data$source$xmlpath[selBrukerXML], selBrukerClass, T)$name )
      dataPaths[[i]]$filepath <- procParams$data$source$xmlpath[selBrukerXML]
      dataPaths[[i]]$brukerclass <- selBrukerClass
    }
    else
    {
      dataPaths[[i]]$name <- basename(procParams$data$source$datapath[i])
      dataPaths[[i]]$filepath <-  procParams$data$source$datapath[i]
    }
  }
  
  #Ask the user for the XML's containing the ROI files
  imgs_names <- unlist(lapply(dataPaths, function(x){ x$name }))
  xmlRoiFiles <- XmlRoiSelectionDialog(imgs_names, init_dir = dirname(procParams$data$source$datapath) )
  if(is.null(xmlRoiFiles))
  {
    #Process aborted by user
    cat("Processing aborted\n")
    return()
  }
  
  for( i in 1:NumOfImages)
  {
    if(xmlRoiFiles$xml_include[i] != "")
    {
      dataPaths[[i]]$xmlroiInclude <- xmlRoiFiles$xml_include[i]
    }
    
    if(xmlRoiFiles$xml_exclude[i] != "")
    {
      dataPaths[[i]]$xmlroiExclude <- xmlRoiFiles$xml_exclude[i]
    }
  }
  rm(NumOfImages)
  
  pixelID_list <- list()
  for( i in 1:NumOfProcRuns)
  {
    cat(paste("Processing execution", i, "of", NumOfProcRuns, "\n"))
    
    #Load each image
    mImg_list <- list()
    for( iimg in 1:NumOfImages2Merge)
    {
      if(procParams$mergedatasets)
      {
        loadImgIndex <- iimg
      }
      else
      {
        loadImgIndex <- i
      }
      
      if( procParams$data$source$type == "xmass" )
      {
        mImg_list[[iimg]] <- rMSI::importBrukerXmassImg( procParams$data$source$datapath, procParams$data$pixelsize,  dataPaths[[loadImgIndex]]$filepath, procParams$data$source$spectrumpath, selected_img = dataPaths[[loadImgIndex]]$brukerclass )
      }
      else
      {
        mImg_list[[iimg]]  <- rMSI::LoadMsiData(dataPaths[[loadImgIndex]]$filepath, ff_overwrite = overwriteRamdisk )  
      }
      
      #Include ID's in include ROI's
      if(!is.null( dataPaths[[loadImgIndex]]$xmlroiInclude))
      {
        #Id's are sorted for better performance
        pixelID_list[[iimg]] <- sort(unique(unlist(lapply(rMSI::ReadBrukerRoiXML(mImg_list[[loadImgIndex]], dataPaths[[loadImgIndex]]$xmlroiInclude), function(x){ x$id }))))
      }
      else
      {
        pixelID_list[[iimg]] <- 0
      }
      
      #Excude ID's in exclude ROI's
      if(!is.null( dataPaths[[loadImgIndex]]$xmlroiExclude))
      {
        #Id's are sorted for better performance
        excudeIds <- sort(unique(unlist(lapply(rMSI::ReadBrukerRoiXML(mImg_list[[loadImgIndex]], dataPaths[[loadImgIndex]]$xmlroiExclude), function(x){ x$id }))))
        
        if( pixelID_list[[iimg]][1] == 0 )
        {
          pixelID_list[[iimg]] <- 1:nrow(mImg_list[[loadImgIndex]]$pos) #Pre-load all Id's in case that no include ROI is specified
        }
        
        excludeIdsLocations <- c()
        for( idExcluded in excudeIds )
        {
          exIdloc <- which(pixelID_list[[iimg]] == idExcluded)
          if(length(exIdloc) > 0)
          {
            excludeIdsLocations <- c(excludeIdsLocations, exIdloc)
          }
        }
        
        if(length(excludeIdsLocations) > 0)
        {
          pixelID_list[[iimg]] <- pixelID_list[[iimg]][-excludeIdsLocations]
        }
      }
    }
    
    # Merge data and resample it if mass axis is different and ID pixel selection if xml are provided
    mImg_list <- MergerMSIDataSets(mImg_list, procParams$data$outpath, pixel_id = pixelID_list ) 
    
    #Independently of the selected norms, set also the one used for data summary export
    if(xmlRoiFiles$summary_export)
    {
      if(xmlRoiFiles$summary_norm == "TIC")
      {
        procParams$spectraNormalization$TIC <- T
      }
      else if(xmlRoiFiles$summary_norm == "MAX")
      {
        procParams$spectraNormalization$MAX <- T
      }
      else if(xmlRoiFiles$summary_norm == "RMS")
      {
        procParams$spectraNormalization$RMS <- T
      }
      else if(xmlRoiFiles$summary_norm == "AcqTic")
      {
        procParams$spectraNormalization$AcqTIC <- T
      }
      procParams$spectraNormalization$enabled <- any( procParams$spectraNormalization$TIC, 
                                                      procParams$spectraNormalization$MAX,
                                                      procParams$spectraNormalization$RMS,
                                                      procParams$spectraNormalization$AcqTIC)
    }
    
    #Process data
    procData <- ProcessImage(img = mImg_list,
                             EnableSmoothing = procParams$smoothing$enabled, SmoothingKernelSize = procParams$smoothing$sgkernsize,
                             EnableAlignment = procParams$alignment$enabled, AlignmentIterations = procParams$alignment$iterations, AlignmentMaxShiftppm = procParams$alignment$maxshift,
                             AlignmentBilinear = procParams$alignment$bilinear, AlignmentOversampling = procParams$alignment$oversampling,
                             AlignmentRefLow = procParams$alignment$reflow, AlignmentRefMid = procParams$alignment$refmid, AlignmentRefHigh = procParams$alignment$refhigh,
                             EnableCalibration = procParams$calibration$enabled, CalibrationPeakWin = procParams$calibration$winsize,
                             EnablePeakPicking = procParams$peakpicking$enabled, SNR = procParams$peakpicking$snr, peakWindow = procParams$peakpicking$winsize, peakUpSampling = procParams$peakpicking$oversample,
                             UseBinning = T, BinTolerance = procParams$peakpicking$bintolerance, BinFilter = procParams$peakpicking$binfilter, BinToleranceUsingPPM = procParams$peakpicking$binUsingPPM,
                             EnableSpectraNormalization = procParams$spectraNormalization$enabled, EnableTICNorm = procParams$spectraNormalization$TIC, EnableRMSNorm = procParams$spectraNormalization$RMS, EnableMAXNorm = procParams$spectraNormalization$MAX, EnableTICAcqNorm = procParams$spectraNormalization$AcqTIC,
                             NumOfThreads = procParams$nthreads, CalSpan = calibrationSpan )

    #Store data summary (before than storing img's because the rMSI objects will be removed)
    if(xmlRoiFiles$summary_export)
    {
      cat("Storing data summary files...\n")
      if(!procParams$peakpicking$enabled)
      {
        procData$peakMat<-NULL
      }
      if(procParams$mergedatasets)
      {
        #Use the full XML file list
        xmlList<-xmlRoiFiles$xml_include
      }
      else
      {
        #Use the current XML
        xmlList<-xmlRoiFiles$xml_include[i]
      }
      
      ExportROIAveragesMultiple(procData$procImg, xmlList, procData$peakMat, procParams$data$outpath, xmlRoiFiles$summary_norm)
    }
    
    #Store peak matrix
    if( procParams$peakpicking$enabled )
    {
      if(procParams$mergedatasets)
      {
        pkMatName <- "mergeddata"
      }
      else
      {
        pkMatName <- sub('\\..*$', '',procData$procImg[[1]]$name) #Remove extensions of image name
      }
      StorePeakMatrix( file.path(procParams$data$outpath,  paste(pkMatName,"-peaks.zip", sep ="")), procData$peakMat)
    }
    
    #Store MS image to a tar file
    if( procParams$smoothing$enabled || procParams$alignment$enabled || procParams$calibration$enabled || procParams$spectraNormalization$enabled || procParams$data$source$type != "tar"  )
    {
      
      for( iSave in 1:length(procData$procImg) )
      {
        imgName <- sub('\\..*$', '',procData$procImg[[iSave]]$name) #Remove extensions of image name
        rMSI::SaveMsiData( procData$procImg[[iSave]], file.path(procParams$data$outpath,  paste(imgName,"-proc.tar", sep ="")))
        
        #Delete data 
        if(deleteRamdisk)
        {
          rMSI::DeleteRamdisk(procData$procImg[[iSave]])
        }
      }
      rm(mImg_list)
      gc()
    }
  }
  
  #Store processing params
  SaveProcessingParams( procParams, 
                        file.path(procParams$data$outpath, "processing_parameters.txt"),
                        xmlRoiFilesInclude = xmlRoiFiles$xml_include, xmlRoiFilesExclude = xmlRoiFiles$xml_exclude,
                        RoiNormalization = xmlRoiFiles$summary_norm)
  
  cat(paste0("All data processing was completed at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S") ,"\n"))
}

#' SaveProcessingParams.
#' 
#' Save all parameters in a list of processing params generated using ImportWizardGui() function.
#' Parameters will be saved in a plain text file.
#'
#' @param procParams a list of parameters.
#' @param filepath a full path where params will be stored
#' @param xmlRoiFilesInclude a vector with the used ROI XML files for ID inclusion, NULL if no ROI was used.
#' @param xmlRoiFilesExclude a vector with the used ROI XML files for ID exclusion, NULL if no ROI was used.
#' @param RoiNormalization a string with the name of normalization used to export the data summary.
#'
SaveProcessingParams <- function( procParams, filepath, xmlRoiFilesInclude = NULL, xmlRoiFilesExclude = NULL, RoiNormalization = NULL)
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
    writeLines(paste("Alignment Bilinear = ", procParams$alignment$bilinear, sep ="" ), con = fObj)
    writeLines(paste("Alignment Oversampling = ", procParams$alignment$oversampling, sep ="" ), con = fObj)
    writeLines(paste("Alignment Ref. Low = ", procParams$alignment$reflow, sep ="" ), con = fObj)
    writeLines(paste("Alignment Ref. Mid = ", procParams$alignment$refmid, sep ="" ), con = fObj)
    writeLines(paste("Alignment Ref. High = ", procParams$alignment$refhigh, sep ="" ), con = fObj)
  }
  
  writeLines(paste("Calibration enabled = ", procParams$calibration$enabled, sep ="" ), con = fObj)
  if(procParams$calibration$enabled)
  {
    #Reserverd for future expansion...
  }

  if(procParams$spectraNormalization$enabled)
  {
    writeLines(paste("Spectra TIC normalization enabled = ", procParams$spectraNormalization$TIC, sep ="" ), con = fObj)
    writeLines(paste("Spectra RMS normalization enabled = ", procParams$spectraNormalization$RMS, sep ="" ), con = fObj)
    writeLines(paste("Spectra MAX normalization enabled = ", procParams$spectraNormalization$MAX, sep ="" ), con = fObj)
    writeLines(paste("Spectra TICAcq normalization enabled = ", procParams$spectraNormalization$AcqTIC, sep ="" ), con = fObj)
  }
  
  writeLines(paste("Peak-picking enabled = ", procParams$peakpicking$enabled, sep ="" ), con = fObj)
  if( procParams$peakpicking$enabled)
  {
    writeLines(paste("Peak-picking SNR threshold = ", procParams$peakpicking$snr, sep ="" ), con = fObj)
    writeLines(paste("Peak-picking detector window = ", procParams$peakpicking$winsize, sep ="" ), con = fObj)
    writeLines(paste("Peak-picking oversampling = ", procParams$peakpicking$oversample, sep ="" ), con = fObj)
    if(procParams$peakpicking$binUsingPPM)
    {
      writeLines("Peak-picking binning tolerance is in ppm", con = fObj)  
    }
    else
    {
      writeLines("Peak-picking binning tolerance is in scans", con = fObj)  
    }
    writeLines(paste("Peak-picking binning tolerance = ", procParams$peakpicking$bintolerance, sep ="" ), con = fObj)
    writeLines(paste("Peak-picking binning filter = ", procParams$peakpicking$binfilter, sep ="" ), con = fObj)
  }
  
  writeLines(paste("Merge datasets = ", procParams$mergedatasets, sep ="" ), con = fObj)
  
  if(!is.null(xmlRoiFilesInclude))
  {
    writeLines("\nUsed ROI Include files:", con = fObj)
    for(i in 1:length(xmlRoiFilesInclude))
    {
      writeLines(paste0("ROI", i, ": ", xmlRoiFilesInclude[i] ), con = fObj)
    }
    writeLines("\n", con = fObj)
  }
  
  if(!is.null(xmlRoiFilesExclude))
  {
    writeLines("\nUsed ROI Exclude files:", con = fObj)
    for(i in 1:length(xmlRoiFilesExclude))
    {
      writeLines(paste0("ROI", i, ": ", xmlRoiFilesExclude[i] ), con = fObj)
    }
    writeLines("\n", con = fObj)
  }
  
  if(!is.null(RoiNormalization))
  {
    writeLines(paste0("Data summary normalization: ", RoiNormalization ), con = fObj)
  }
  
  close(fObj)
}

#' MergerMSIDataSets.
#' 
#' Merges various rMSI objects in order to process all of them in the same run.
#' For example, this is usefull to align data form diferent experiments together.
#' In case that images don't share the same mass axis, all of them will be resampled and stored in a new ramdisk.
#' If pixel_id is provided, a new ramdisk for each image will be created to filter out non specified pixel ID's.
#' Discarding pixels of a dataset may result interesting to avoid artifacts in alignment and peak binning if some
#' regions are not well-correlated to the rest of tissue.
#'
#' @param img_list a list of images to be merged.
#' @param ramdisk_path a path where resampled data ramdisk will be stored.
#' @param pixel_id a list containing a vector of ID's to retain for each img. If some img have to use all ID's 0 my be supplied. 
#'
#' @return a list with the merged images.
#'
MergerMSIDataSets <- function( img_list, ramdisk_path, pixel_id = NULL )
{
  # 1- If the img_list contains only one element, then is nothing to merge, just return the list itself
  if( length(img_list) == 1 && is.null(pixel_id))
  {
    return( img_list )
  }
  if( length(img_list) == 1 && pixel_id[[1]] == 0)
  {
    return( img_list )
  }
  
  # 2- Check if the mass axis is the same
  common_mass <- img_list[[1]]$mass
  identicalMassAxis <- TRUE
  if( length(img_list) > 1)
  {
    for( i in 2:length(img_list))
    {
      identicalMassAxis <- identicalMassAxis & identical(common_mass, img_list[[i]]$mass)
    }
  }
  
  if(identicalMassAxis && is.null(pixel_id))
  {
    #No need to merge, all data sets can be used together because they share the same mass axis
    return( img_list )
  }
  
  if(identicalMassAxis &&  all( unlist(lapply(pixel_id, function(x){ x == 0}))))
  {
    #No need to merge, all data sets can be used together because they share the same mass axis and no ROI filtering used
    return( img_list ) 
  }
  
  # 3- Calculate the new common mass axis
  common_mass <- img_list[[1]]$mass
  if(!identicalMassAxis)
  {
    for( i in 2:length(img_list))
    {
      common_mass <- rMSI::MergeMassAxis(common_mass, img_list[[i]]$mass)
    }
  }
  
  # 4- Apply the resampling and/or pixel id selection and create the new ramdisk 
  for( i in 1:length(img_list))
  {
    cat(paste0("Merging dataset ", i, "/", length(img_list), "\n"))
    
    # Get ID's for current image
    if(is.null(pixel_id))
    {
      imgSelIDs <- 1:nrow(img_list[[i]]$pos)
    }
    else if ( pixel_id[[i]][1] == 0 )
    {
      imgSelIDs <- 1:nrow(img_list[[i]]$pos)
    }
    else
    {
      imgSelIDs <-  pixel_id[[i]]
    }

    data_path <- file.path(ramdisk_path, paste0(img_list[[i]]$name, "_resam"))
    imgAux <- rMSI::CreateSubDataset(img_list[[i]], imgSelIDs, data_path, common_mass)
    rMSI::DeleteRamdisk(img_list[[i]])
    img_list[[i]] <- imgAux
  }
  
  return(img_list)
}
