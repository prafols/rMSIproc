#########################################################################
#     rMSIproc - R package for MSI data processing
#     Copyright (C) 2017 Pere Rafols Soler
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

#' ExportROIAverages.
#' 
#' Exports the average spectrum of a group of pixels in a ROI as an ASCII text file.
#' ROI's are defined using Bruker XML files. A single spectrum will be exported for each ROI found in XML file.
#'
#' @param roi_xml_file a Bruker region XML file.
#' @param img an rMSI object.
#' @param out_path  disk path where the results will be stored.
#' @param normalization a text string with the name of desired normalization to apply at output files.
#'
#' @return a list with the data set UUID and the pixels ID's exported for each ROI.
#'
ExportROIAverages <- function( roi_xml_file, img, out_path = getwd(), normalization = NULL)
{
  #Read the normalization factor
  if( !is.null(normalization) )
  {
    if(is.null(img$normalizations))
    {
      stop("Error: No normalizations found in rMSI object.\n")
    }
    if(is.null(img$normalizations[[normalization]]))
    {
      stop("Error: The specified normalzation is not available in rMSI object.\n")
    }
    norm_factors <- img$normalizations[[normalization]]
  }
  else
  {
    #Normalization disabled
    norm_factors <- rep(1, nrow(img$pos))
  }
  
  #Prepara the output paths
  asciiSpcPath <- file.path(out_path, "roi_spectra")
  dir.create(asciiSpcPath, showWarnings = F, recursive = T)
  
  # Get the rois listing IDs  
  id_list <- rMSI::ReadBrukerRoiXML( img, roi_xml_file )
  
  # Iterate for each roi
  for( i in 1:length(id_list))
  {
    #Calculate average spectrum
    cat(paste0("Exporting average spectrum of roi ", i, "/", length(id_list), "\n"))
    cubes <- rMSI::getCubeRowFromIds(img, id_list[[i]]$id) 
    averageSpc <- rep(0, length(img$mass))
    for( ic in 1:length(cubes))
    {
      dm <- img$data[[cubes[[ic]]$cube]][cubes[[ic]]$row, ]
      dm <- dm/norm_factors[cubes[[ic]]$id] #Apply normalization
      averageSpc <- averageSpc + colSums(dm)
    }
    averageSpc <- averageSpc / length( id_list[[i]]$id) #Div to total number of pixels in ROI
    
    #Store the average spectrum to a ASCII file
    write.table(matrix(data = c(img$mass, averageSpc),  ncol = 2, nrow = length(averageSpc)), 
                file = file.path(asciiSpcPath, paste0(img$name, "_", id_list[[i]]$name, ".txt")), 
                append = F, row.names = F, col.names = F, sep = "\t", eol = "\r\n")
  }
  
  return( list( img_uuid = img$uuid, idlist = id_list ))
}

#' ExportROIPeaks.
#' 
#' Export a peak matrix summary as CSV files.
#' An average and standard deviation of each ROI is calculated for intensity, area and SNR peak matrices.
#' The results are stored as tables in CSV files.
#'
#' @param pkmat an rMSIproc peak matrix object.
#' @param idlist a list of pixel ID's in each ROI using the format returned by ExportROIAverages.
#' @param out_path disk path where the results will be stored.
#' @param normalization a text string with the name of desired normalization to apply at output files.
#'
ExportROIPeaks <- function(pkmat, idlist, out_path = getwd(), normalization = NULL)
{
  #Read the normalization factor
  if( !is.null(normalization) )
  {
    if(is.null(pkmat$normalizations))
    {
      stop("Error: No normalizations found in rMSIproc peak matrix.\n")
    }
    if(is.null(pkmat$normalizations[[normalization]]))
    {
      stop("Error: The specified normalzation is not available in rMSIpeak matrix\n")
    }
    norm_factors <- pkmat$normalizations[[normalization]]
  }
  else
  {
    #Normalization disabled
    norm_factors <- rep(1, nrow(pkmat$intensity))
  }
  
  #Count total number of rows in final data.frame to track a progress bar and prepare the data.frame
  roiNames <- unlist(lapply(idlist, function(x){ unlist(lapply(x$idlist, function(y){ y$name} )) }))
  numOfRois <- length(roiNames)
  roiNames <- c("ROI", roiNames) #The first row of the data.frame is the column names
  emtyAuxMat <- matrix(ncol = length(pkmat$mass), nrow = numOfRois + 1) #+1 because the first row will be used for mass axis in numeric format
  colnames(emtyAuxMat) <- sprintf("mz%.6f", pkmat$mass)
  emtyAuxMat[1, ] <- pkmat$mass
  dfintAvg <- data.frame( ROI = roiNames , emtyAuxMat)
  dfAreaAvg <- data.frame( ROI = roiNames , emtyAuxMat)
  dfSNRAvg <- data.frame( ROI = roiNames , emtyAuxMat)
  dfintSD <- data.frame( ROI = roiNames , emtyAuxMat)
  dfAreaSD <- data.frame( ROI = roiNames , emtyAuxMat)
  dfSNRSD <- data.frame( ROI = roiNames , emtyAuxMat)
  rm(emtyAuxMat)
  
  pb <- txtProgressBar(min=0, max=numOfRois, style = 3)
  currentROi <- 0
  for( img in 1:length(idlist))
  {
    #Select the appropiate pkMat corresponding to current image
    selpkMat <- which(pkmat$uuid == idlist[[img]]$img_uuid)
    if(length(selpkMat) == 0)
    {
      stop("Error: The name of rMSI images does not match any object in rMSIproc peak matrix.\n")
    }
    
    for( iroi in 1:length(idlist[[img]]$idlist))
    {
      currentROi <- currentROi + 1
      setTxtProgressBar(pb, currentROi)
      
      roi_rows <- rMSIproc::getPeakMatrixRowsFromImgIds(pkmat, selpkMat, idlist[[img]]$idlist[[iroi]]$id)
      dmInt <- pkmat$intensity[roi_rows, ] / norm_factors[roi_rows]
      dmArea <- pkmat$area[roi_rows, ] / norm_factors[roi_rows]
      dmSNR <- pkmat$SNR[roi_rows, ]
      
      #Average rows
      dfintAvg[currentROi + 1, 2:(1+length(pkmat$mass))] <- apply(dmInt, 2, mean)
      dfAreaAvg[currentROi + 1, 2:(1+length(pkmat$mass))] <- apply(dmArea, 2, mean)
      dfSNRAvg[currentROi + 1, 2:(1+length(pkmat$mass))] <- apply(dmSNR, 2, mean)
      
      #Standard dev. rows
      dfintSD[currentROi + 1, 2:(1+length(pkmat$mass))] <- apply(dmInt, 2, sd)
      dfAreaSD[currentROi + 1, 2:(1+length(pkmat$mass))] <- apply(dmArea, 2, sd)
      dfSNRSD[currentROi + 1, 2:(1+length(pkmat$mass))] <- apply(dmSNR, 2, sd)
    }
  }
  close(pb)
  
  #Store data to CSV
  dir.create(file.path(out_path, "roi_peaks"), showWarnings = F, recursive = T)
  
  if(length(pkmat$names)>1)
  {
    baseFname <- "merged" 
  }
  else
  {
    baseFname <- pkmat$names[1]
  }
  
  write.table(dfintAvg, file = file.path(out_path, "roi_peaks", paste0(baseFname, "_intensity_average.csv")), eol = "\r\n", row.names = F, col.names = F, dec = ".", sep = ";" )
  write.table(dfAreaAvg, file = file.path(out_path, "roi_peaks", paste0(baseFname, "_area_average.csv")), eol = "\r\n", row.names = F, col.names = F, dec = ".", sep = ";")
  write.table(dfSNRAvg, file = file.path(out_path, "roi_peaks", paste0(baseFname, "_snr_average.csv")), eol = "\r\n", row.names = F, col.names = F, dec = ".", sep = ";")
  write.table(dfintSD, file = file.path(out_path, "roi_peaks", paste0(baseFname, "_intensity_sdev.csv")), eol = "\r\n", row.names = F, col.names = F, dec = ".", sep = ";")
  write.table(dfAreaSD, file = file.path(out_path, "roi_peaks", paste0(baseFname, "_area_sdev.csv")), eol = "\r\n", row.names = F, col.names = F, dec = ".", sep = ";")
  write.table(dfSNRSD, file = file.path(out_path, "roi_peaks", paste0(baseFname, "_snr_sdev.csv")), eol = "\r\n", row.names = F, col.names = F, dec = ".", sep = ";")
}

#' ExportROIAveragesMultiple.
#' 
#' Exports the ROI average spectrum for a data set of multiple images. If pkMat != NULL the CSV summary peak
#' matrices will be exported as well.
#'
#' @param img_list a lis of rMSI objects.
#' @param xml_list a vector of XML ROI files in the same order as img_list.
#' @param pkMat an rMSIproc peak matrix object.
#' @param out_path disk path where the results will be stored.
#' @param normalization a text string with the name of desired normalization to apply at output files.
#'
ExportROIAveragesMultiple <- function( img_list, xml_list, pkMat = NULL, out_path = getwd(), normalization = NULL)
{
  if(length(img_list) != length(xml_list))
  {
    stop("Error: img list and xml list must have the same length.\n")
  }
  
  idLists <- list()
  for( i in 1:length(img_list))
  {
    idLists[[i]]<-ExportROIAverages(xml_list[i], img_list[[i]], out_path , normalization) 
  }
  
  if(!is.null(pkMat))
  {
    ExportROIPeaks(pkMat, idLists, out_path, normalization)
  }
}
