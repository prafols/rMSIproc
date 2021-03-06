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

#Wrapper to move ff data objects in rMSI format between R and C

#' printrMSIdataInfo.
#'
#' Prints HDD storing information of rMSI data object using the C++ backend.
#' 
#' @param img and rMSI object image.
#' @export
#'
printrMSIdataInfo <- function( img )
{
  dataInf <- getrMSIdataInfo(img)
  PrintrMSIObjectInfo( dataInf$filenames, dataInf$masschannels, dataInf$nrows, dataInf$datatype )  
}

#' loadDataCube.
#' 
#' Loads an rMSI data cube (aka. ramdisk file) to an R matrix using C++ backend.
#' This method is just for testing convinience and should not be used because it copies a matrix indexed by rows to a
#' matrix indexed by columns so it performs very slowly. Instead of using this method use rMSI::loadImgCunckFromCube.
#'
#' @param img  and rMSI object image.
#' @param cubeSel the cube to load R index.
#'
#' @return an R matrix with spectra in rows.
#' @export
loadDataCube <- function( img, cubeSel )
{
  dataInf <- getrMSIdataInfo(img)
  return(LoadrMSIDataCube( dataInf$filenames, dataInf$masschannels, dataInf$nrows, dataInf$datatype, cubeSel - 1))
}


#' getrMSIdataInfo.
#' 
#' Obtains all storing information of an rMSI object and returns it as a list.
#'
#' @param img an rMSI data object.
#'
#' @return a list containing all storing information.
#' 
getrMSIdataInfo <- function( img )
{
  ffData <- list()
  ffData$filenames <- unlist(lapply(img$data, function(x){ path.expand(attr(attr(x, "physical"), "filename")) }))
  ffData$masschannels <- length(img$mass)
  ffData$nrows <- unlist(lapply(img$data, function(x) { attr(attr(x, "virtual"), "Dim")[1]}))
  ffData$datatype <- attr(attr(img$data[[1]], "physical"), "vmode") 
  return(ffData)
}

#' getrMSIdataInfoMultipleDataSets.
#' 
#' Obtains all storing information of various rMSI objects and returns them unified in a list.
#' The supplied rMSI objects must have the same number of mass channels and the same data type.
#' Otherwise an error will be raised.
#' The returned list contains an extra file "datasets" to indecate at which dataset belongs each ramdisk.
#'
#' @param imgs_list a list of rMSI objects.
#'
#' @return a list containing all storing information unified. 
#' 
getrMSIdataInfoMultipleDataSets <- function( imgs_list )
{
  imgInfo <- getrMSIdataInfo(imgs_list[[1]])
  imgInfo$dataset <- rep(1, length( imgInfo$filenames ))
  
  if( length(imgs_list) > 1)
  {
    for( i in 2:length(imgs_list))
    {
      aux <- getrMSIdataInfo(imgs_list[[i]])
      
      if( imgInfo$masschannels != aux$masschannels )
      {
        stop("Error: Number of mass channels of images are different\n")
      }
      
      if( !identical(imgs_list[[1]]$mass, imgs_list[[1]]$mass) )
      {
        stop("Error: Mass axis of images are different\n")
      }
      
      if( imgInfo$datatype != aux$datatype)
      {
        stop("Error: Data type of images are different\n")
      }
      
      imgInfo$filenames <- c(imgInfo$filenames, aux$filenames)
      imgInfo$nrows <- c( imgInfo$nrows,  aux$nrows)
      imgInfo$dataset <- c( imgInfo$dataset, rep(i, length( aux$filenames )))
    }
  }
  
  return(imgInfo)
}
