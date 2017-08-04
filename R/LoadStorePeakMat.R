#' LoadPeakMatrix.
#' 
#' Loads a binned peaks matrix from HDD.
#'
#' @param data_path full path to zip file where data is stored.
#'
#' @return an R List containing intensity, SNR and area matrices, mass axis vector and if available the normalizations data.frame.
#'
LoadPeakMatrix <- function( data_path )
{
  dir.create(file.path(dirname(data_path), "tmp"), recursive = T)
  cat("Unzipping data...\n")
  unzip(data_path, overwrite = T, exdir = file.path(dirname(data_path), "tmp"))
  
  if(dir.exists(file.path(dirname(data_path), "tmp")))
  {
    ldata <-  LoadPeakMatrixC( file.path(dirname(data_path), "tmp") )
  }
  else
  {
    stop("Error: the specified data dire does not exist.\n")
  }
  
  unlink(file.path(dirname(data_path), "tmp"), recursive = T)
  class(ldata) <- "rMSIprocPeakMatrix"
  
  colnames(ldata$pos) <- c("x", "y")
  if(!is.null(ldata$posMotors))
  {
    colnames(ldata$posMotors) <- c("x", "y")
  }
  
  cat("Done\n")
  return(ldata)
}

#' StorePeakMatrix.
#'
#' Stores a binned peaks matrix to HDD.
#' Data is stored zip compressed, so it is recomeneded to specify the name with .zip extension.
#'
#' @param data_path full path including filename where data must be stored.
#' @param data a List containing intensity, SNR and area matrices, the mass axis vector and a data.frame containing in each variable a normalization vector. 
#'
StorePeakMatrix <- function( data_path, data )
{
  if( class(data) != "rMSIprocPeakMatrix")
  {
    stop("The provided peak matrix is not in rMSIprocPeakMatrix format\n")
  }
  
  dir.create(file.path(dirname(data_path), "tmp"), recursive = T)
  StorePeakMatrixC( file.path(dirname(data_path), "tmp"), data )
  cat("Adding data to zip archive...\n")
  old_wd <- getwd()
  setwd(file.path(dirname(data_path), "tmp"))
  zip(data_path, files = dir())
  setwd(old_wd)
  unlink(file.path(dirname(data_path), "tmp"), recursive = T)
  cat("Done\n")
}

#' buildImgIdVectorFromPeakMatrix.
#' 
#' Builds a integer vector containing all rMSI objects ID's accroding the peak matrix row order.
#' The resulting ID vector can be used to locate a spectrum ID in a peak matrix of multiple data sets.
#'
#' @param pkMat an rMSIproc peak matrix object.
#'
#' @return an integer vector with all ID's of each image in the pkMat object.
#'
buildImgIdVectorFromPeakMatrix <- function(pkMat)
{
  imgIDs <- rep(NA, sum(pkMat$numPixels))
  istart <- 1
  for( i in 1:length(pkMat$numPixels))
  {
    istop <- istart + pkMat$numPixels[i] - 1
    imgIDs[istart:istop] <- 1:pkMat$numPixels[i]
    istart <- istop + 1
  }
  return(as.integer(imgIDs))
}

#' getImgIdsFromPeakMatrixRows.
#' 
#' Calculate the rMSI object ID's corresponding to a subset of row of a rMSIproc peak matrix.
#'
#' @param pkMat an rMSIproc peak matrix object.
#' @param rows a vector of peak matrix rows.
#'
#' @return a list with the rMSI object ID's and image that correpond the specified rows.
#'
getImgIdsFromPeakMatrixRows <- function(pkMat, rows)
{
  #Sort rows assending and remove duplicates for fater performance
  rows <- sort(unique(rows), decreasing = F)
  
  id_lst <- list()
  allIds <- buildImgIdVectorFromPeakMatrix(pkMat)
  iRowStart <- 1
  for( i in 1:length(pkMat$names))
  {
    iRowStop <- iRowStart + pkMat$numPixels[i] - 1
    iSubRows<- which(rows %in% iRowStart:iRowStop )
    if( length(iSubRows) > 0)
    {
      id_lst[[i]] <- list( name = pkMat$names[i], id = allIds[rows[iSubRows]] , pkMatRow = rows[iSubRows])
    }
    else
    {
      id_lst[[i]] <- list( name = pkMat$names[i], id = c() , pkMatRow = c())
    }
    iRowStart <- iRowStop + 1
  }
  return(id_lst)
}

#' getPeakMatrixRowsFromImgIds.
#' 
#' Obtains a vector of rMSIproc peak matrix rows corresponding to a vector of rMSI obj ID's.
#'
#' @param pkMat an rMSIproc peak matrix object.
#' @param img_num the number of the data set object to select in pkMat.
#' @param ids a vector of rMSI obj ID's.
#'
#' @return a vector containing the rows of peak matrix that correspond to the selected rMSI object ID's.
#'
getPeakMatrixRowsFromImgIds <- function( pkMat, img_num, ids )
{
  if(length(pkMat$numPixels) < img_num)
  {
    stop(paste("Error: pkMat does not containg", img_num, "images.\n"))
  }
  
  if(img_num > 1)
  {
    startRow <- sum(pkMat$numPixels[1:(img_num-1)]) + 1 
  }
  else
  {
    startRow <- 1
  }
  stopRow <- startRow + pkMat$numPixels[img_num] - 1
  selRows <- ids + startRow - 1
  
  if(max(selRows) > stopRow || min(selRows) < startRow)
  {
    stop("Error: Id selection is out of range.\n")
  }
  
  return(selRows)
}
