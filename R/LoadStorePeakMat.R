#' LoadPeakMatrix.
#' 
#' Loads a binned peaks matrix from HDD.
#'
#' @param data_path full path to zip file where data is stored.
#'
#' @return an R List containing intensity, SNR and area matrices and mass axis vector.
#' @export
#'
LoadPeakMatrix <- function( data_path)
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
  
  cat("Done\n")
  return(ldata)
}

#' StorePeakMatrix.
#'
#' Stores a binned peaks matrix to HDD.
#' Data is stored zip compressed, so it is recomeneded to specify the name with .zip extension.
#'
#' @param data_path full path including filename where data must be stored.
#' @param data a List containing intensity, SNR and area matrices and mass axis vector.
#'
#' @export
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