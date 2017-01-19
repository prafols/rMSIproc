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
#' The reference spectrum is the spectrum in the dataset with the bes correlation to average spectrum.
#'
#' @param img MSI image to calculate internal reference spectrum.
#'
#' @return an intensity vector corresponding to the reference spectrum.
#' @export
#'
InternalReferenceSpectrum <- function(img)
{
  cat("Calculating internal reference spectrum...\n")
  #pb <- txtProgressBar(min = 0,  max = nrow(img$pos), style = 3) #TODO disabled for best testing
  id <- 1
  maxCor <- 0
  maxId <- 0
  for( i in 1:length(img$data))
  {
    cat(paste("Cube:",i, "Row:", j, "\n"))
    dc <- rMSI::loadImgCunckFromCube(img, i)
    for( j in 1:nrow(dc))
    {
      pxCor <- cor(img$mean, dc[j, ] )
      if( pxCor > maxCor )
      {
        maxCor <- pxCor
        maxId <- id
      }
      #setTxtProgressBar(pb, id)  #TODO disabled for best testing
      id <- id + 1
    }
  }
  #close(pb)  #TODO disabled for best testing
  return(rMSI::loadImgCunckFromIds(img, maxId)[1,])
}
