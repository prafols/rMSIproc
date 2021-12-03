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

#' InternalReferenceSpectrum.
#' 
#' Calculates the dataset reference spectrum to use in label free alignment.
#' The reference spectrum is the spectrum in the dataset with the best correlation to average spectrum and with high TIC.
#'
#' @param img MSI image to calculate internal reference spectrum.
#' @param reference the spectrum to use as reference for correlations.
#'
#' @return a list with the intensity vector corresponding to the reference spectrum, the score and the pixel ID selected as reference.
#'
InternalReferenceSpectrum <- function(img, reference = img$mean)
{
  pb <- txtProgressBar(min = 0,  max = nrow(img$pos), style = 3)
  id <- 1
  maxScore <- 0
  maxId <- 0
  selID <- NA
  TICref <- sum(reference)
  for( i in 1:length(img$data))
  {
    dc <- rMSI:::ffVector2Matrix(rMSI::loadImgChunkFromCube(img, i))
    for( j in 1:nrow(dc))
    {
      if(var(dc[j, ]) > 0)
      {
        pxCor <- cor(reference, dc[j, ] )
        ticScr <- sum(dc[j, ])/TICref
        score <- pxCor*ticScr
        if( score > maxScore )
        {
          maxScore <- score
          maxId <- id
          selID <- id
        }
      }
      setTxtProgressBar(pb, id)
      id <- id + 1
    }
  }
  
  close(pb)
  return( list(spectrum = rMSI::loadImgChunkFromIds(img, maxId)[1,], score = maxScore, ID = selID ))
}

#' InternalReferenceSpectrumMultipleDatasets.
#' 
#' Calculates the dataset reference spectrum to use in label free alignment.
#' The reference spectrum is the spectrum in the dataset with the best correlation to average spectrum and with high TIC.
#'
#' @param img_list a list of various rMSI objects.
#'
#' @return a list with the intensity vector corresponding to the reference spectrum, the score, the image index which contains the refernce spectrum and the pixel ID selected as reference.
#' 
InternalReferenceSpectrumMultipleDatasets <- function(img_list)
{
  #Set the global reference as the average of the first dataset
  globalAverage <- img_list[[1]]$mean 
  
  #Calculate correlations
  bestRefs <- list()
  for( i in 1:length(img_list) )
  {
    cat(paste0("Calculating internal reference spectrum ", i, "/",length(img_list),"...\n"))
    bestRefs[[i]] <- InternalReferenceSpectrum(img_list[[i]], globalAverage)
  }
  
  #Select the best correlation in all datasets
  bestID <- which.max(unlist(lapply(bestRefs, function(x){ return(x$score) })))

  return(  list(spectrum = bestRefs[[bestID]]$spectrum, score = bestRefs[[bestID]]$score, imgIndex = bestID, ID =  bestRefs[[bestID]]$ID ) )
}
