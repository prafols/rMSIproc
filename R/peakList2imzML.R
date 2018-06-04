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


#' export_imzMLpeakList.
#' 
#' Export an rMSIproc peak list as an imzML processed dataset.
#'
#' @param peakList an rMSIproc peak list object.
#' @param posMatrix an rMSI pos matrix with the pixel coordinates in the image.
#' @param pixel_size_um the image pixel size in microns.
#' @param filename complete path where the imzML and ibd files will be stored (the .imzML extensions must be omited).
#' @param normalizations a numeric vector contaning the normalization value for each pixel (1 by default).
#'
#' @export
#'
export_imzMLpeakList <- function(peakList, posMatrix, pixel_size_um, filename, normalization = rep(1,length(peakList)))
{
  if(length(peakList) != nrow(posMatrix))
  {
    stop("peakList and posMatrix must be the same length.\n")
  }
  if(length(peakList) != length(normalization))
  {
    stop("peakList and posMatrix must be the same normalization.\n")
  }
  if(length( which(normalization == 0)) > 0)
  {
    stop("Invalid normalization vector, some elements are zero.\n")   
  }
  
  dir.create(dirname(path.expand(filename)), showWarnings = F, recursive = T)
  baseFileName <- basename(filename)
  
  #1- Create the offset matrix
  offMat <- list()
  offMat$UUID <- rMSI:::uuid_timebased()
  offMat$continuous_mode <- FALSE
  offMat$compression_mz <- FALSE
  offMat$compression_int <- FALSE
  offMat$mz_dataType <- "double"
  offMat$int_dataType <- "double"
  offMat$pixel_size_um <- pixel_size_um
  
  offMat$run_data <- data.frame( x =  posMatrix[,"x"],
                                 y =  posMatrix[,"y"],
                                 mzLength = unlist(lapply(peakList, function(x){length(x$mass)})),
                                 mzOffset = rep(0,nrow(posMatrix)), #Offsets are calculated in the next for loop
                                 intLength = unlist(lapply(peakList, function(x){length(x$intensity)})),
                                 intOffset =  rep(0,nrow(posMatrix)) #Offsets are calculated in the next for loop
  )
  offMat$run_data$mzOffset[1] <- 16
  offMat$run_data$intOffset[1] <- offMat$run_data$mzOffset[1] + 8*offMat$run_data$mzLength[1]
  for( i in 2:nrow(offMat$run_data))
  {
    offMat$run_data$mzOffset[i] <- offMat$run_data$intOffset[i-1] + 8*offMat$run_data$intLength[i-1]
    offMat$run_data$intOffset[i] <- offMat$run_data$mzOffset[i] + 8*offMat$run_data$mzLength[i]
  }
  
  #1- Store the binary data
  cat("Writing the .ibd file(binary)...\n")
  pb <- txtProgressBar(min = 0, max = nrow(offMat$run_data), initial = 0, style = 3)
  ibd_fname <- file.path(dirname(path.expand(filename)), paste0(baseFileName, ".ibd"))
  ibd_conn <- file(ibd_fname, "w+b")
  intUUID <- strtoi(substring(offMat$UUID, seq(1,nchar(offMat$UUID),2), seq(2,nchar(offMat$UUID),2)), base = 16)
  writeBin(intUUID, ibd_conn, size = 1, endian="little") #Write the UUID
  for( i in 1:nrow(offMat$run_data))
  {
    setTxtProgressBar(pb, i)
    writeBin(peakList[[i]]$mass, ibd_conn, size = 8, endian="little", useBytes = T)
    writeBin(peakList[[i]]$intensity/normalization[i], ibd_conn, size = 8, endian="little", useBytes = T)
  }
  close(pb)
  close(ibd_conn)
  
  cat("Calculating MD5 checksum...\n")
  offMat$SHA <- ""
  offMat$MD5 <- toupper(digest::digest( ibd_fname, algo = "md5", file = T))
  
  #2- Store the xml part
  cat("Writing the .imzML file(XML)...\n")
  imzML_fname <- file.path(dirname(path.expand(filename)), paste0(baseFileName, ".imzML"))
  bResult <- rMSI:::CimzMLStore( imzML_fname ,offMat)
  if(bResult)
  {
    cat("imzML exported successfully.\n")
  }
  else
  {
    cat("imzML exported with errors.\n")
  }
}


#' import_imzMLpeakList.
#' 
#' import a processed imzML file to a rMSIproc peak list object.
#'
#' @param imzML_File full path to .imzML file.
#' @param ibd_File path to the binary file (default the same as imzML file but with .ibd extension)
#'
#' @return a list containing: the peakList, the rMSI formated position matrix and the pixel size.
#' @export
#'
import_imzMLpeakList <- function(imzML_File, ibd_File =  paste(sub("\\.[^.]*$", "", imzML_File), ".ibd", sep = "" ))
{
  xmlRes <- rMSI:::CimzMLParse(path.expand(imzML_File))
  if( !is.null(xmlRes$Error))
  {
    stop(paste0(xmlRes$Error, "\n"))
  }
  
  #Check if data is in processed mode
  if(xmlRes$continuous_mode)
  {
    stop("imzML data is in continuous mode. This mode is not supported for peak list data. Use rMSI to load continuous data instead.\n")
  }

  #Checksums
  if( xmlRes$SHA != "" )
  {
    cat("\nChecking binary data checksum using SHA-1 key... ")
    res <- toupper(digest::digest( ibd_File, algo = "sha1", file = T))
    if( res == xmlRes$SHA )
    {
      cat("OK\n")
    }
    else
    {
      cat(paste("NOK\nChecksums don't match\nXML key:", xmlRes$SHA, "\nBinary file key:", res,"\n"))
      cat("WARNING: MS data my be corrupt!\n")
    }
  }
  if( xmlRes$MD5 != "")
  {
    cat("Checking binary data checksum using MD5 key... ")
    res <- toupper(digest::digest( ibd_File, algo = "md5", file = T))
    if( res == xmlRes$MD5 )
    {
      cat("OK\n")
    }
    else
    {
      cat(paste("NOK\nChecksums don't match\nXML key:", xmlRes$MD5, "\nBinary file key:", res,"\n"))
      cat("WARNING: MS data my be corrupt!\n")
    }
  }
  
  #Create a connection to read binary file
  bincon <- file(description = ibd_File, open = "rb")
  
  #Test the UUID in binary file (the first 16 bytes are always UUID (in XML file are in hex codes))
  binUUID <- paste(sprintf("%.2X", readBin(bincon, integer(), 16, size = 1, signed = F)), collapse = "")
  if(binUUID != xmlRes$UUID)
  {
    close(bincon)
    stop("ERROR: UUID in imzML file does not match UUID in ibd file\n")
  }
  
  #Map datatypes
  sizeInBytesFromDataType <- function(str_dataType)
  {
    if(str_dataType == "int" || str_dataType == "float")
    {
      return(4)
    }
    if(str_dataType == "long" || str_dataType == "double")
    {
      return(8)
    }
  }
  if(xmlRes$mz_dataType == "int" || xmlRes$mz_dataType == "long")
  {
    readDataTypeMz <- integer()
  }
  if(xmlRes$mz_dataType == "float" || xmlRes$mz_dataType == "double")
  {
    readDataTypeMz <- numeric()
  }
  bytes2ReadMz <- sizeInBytesFromDataType(xmlRes$mz_dataType)
  if(xmlRes$int_dataType == "int" || xmlRes$int_dataType == "long")
  {
    readDataTypeInt <- integer()
  }
  if(xmlRes$int_dataType == "float" || xmlRes$int_dataType == "double")
  {
    readDataTypeInt <- numeric()
  }
  bytes2ReadInt <- sizeInBytesFromDataType(xmlRes$int_dataType)
  
  #Read bin data
  cat("Reading the .ibd file(binary)...\n")
  pb <- txtProgressBar(min = 0, max = nrow(xmlRes$run_data), initial = 0, style = 3)
  peakLst <- list()
  for( i in 1:nrow(xmlRes$run_data))
  {
    setTxtProgressBar(pb, i)
    peakLst[[i]] <- list()
    seek(bincon, rw = "read", where = xmlRes$run_data[i, "mzOffset"] )
    peakLst[[i]]$mass <- readBin(bincon, readDataTypeMz, xmlRes$run_data[i, "mzLength"], size = bytes2ReadMz, signed = T)
    seek(bincon, rw = "read", where = xmlRes$run_data[i, "intOffset"] )
    peakLst[[i]]$intensity <- readBin(bincon, readDataTypeInt, xmlRes$run_data[i, "intLength"], size = bytes2ReadInt, signed = T)
  }
  close(pb)
  close(bincon)
  
  #The position matrix
  RunCoords <- matrix( c( xmlRes$run_data$x,  xmlRes$run_data$y ), nrow = nrow( xmlRes$run_data), ncol = 2, byrow = F)
  colnames(RunCoords) <- c("x", "y")
  
  return(list(peakList = peakLst, pos = RunCoords, pixel_size_um = xmlRes$pixel_size_um ))
}
