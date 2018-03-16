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

#' OffsetRasterPosMats.
#' 
#' creates an unified pos matrix to plot properly various datasets in a single image.
#'
#' @param posMatrix the pos matrix of an rMSIproc peak matrix object.
#' @param numPixels a vector with the number of pixels of each image.
#' @param horizontal if must be laydout horizontally or vertically.
#' @param margin number in pixels to separe images.
#' 
#' @return a pos Matrix with proper offsets to allow multi dataset visualization
#'
OffsetRasterPosMats <- function(posMatrix, numPixels, horizontal = T, margin = 0)
{
  offsetPosMat <-  posMatrix
  istart <- 1
  offset <- 0
  for( i in 1:length(numPixels))
  {
    istop <- istart + numPixels[i] - 1
    if(horizontal)
    {
      offsetPosMat[istart:istop, "x"] <- offsetPosMat[istart:istop, "x"] + offset + margin
      offset <- max(offsetPosMat[istart:istop, "x"])
    }
    else
    {
      offsetPosMat[istart:istop, "y"] <- offsetPosMat[istart:istop, "y"] + offset + margin
      offset <- max(offsetPosMat[istart:istop, "y"])
    }
    istart <- istop + 1
  }
  return(offsetPosMat)
}

#' plotPeakImage.
#' 
#' plot the ion image map of a given mass or column of a rMSIproc peak matrix object.
#' If the peak matrix contains data from various datasets the images will be layout horizontally or vertically.
#' At leas mz or column must be specified.
#'
#' @param peakMatrix the peak matrix in an rMSIproc object.
#' @param mz the peak mass to plot, the nearest peak mass will be ploted.
#' @param column the column of the peak matrix to plot.
#' @param offset_horizontal if true images will be laydout horizontally.
#' @param matrix the name of the peak matrix to plot.
#' @param normalization the name of normalization to use.
#' @param rotate the rotation in degrees.
#' @param pixel_size_um the pixel resolution in um.
#'
plotPeakImage <- function(peakMatrix, mz = NULL, column = NULL, offset_horizontal = T, matrix = "intensity", normalization = NULL, rotate = 0, pixel_size_um = 100)
{
  if(is.null(peakMatrix[[matrix]]))
  {
    stop("The selected peak matrix does not exist.\n")
  }
  
  if(is.null(mz) && is.null(column))
  {
    stop("At least mass or a column must be specified")
  }
  
  if( !is.null(mz))
  {
    column <- which.min(abs(peakMatrix$mass - mz))
  }
  
  if( !is.null(normalization))
  {
    if( is.null(peakMatrix$normalizations))
    {
      stop("No normalizations avaliable for the supplied peak matrix")
    }
    
    normVals <- peakMatrix$normalizations[[normalization]]
    if(is.null(normVals))
    {
      stop(paste0("No normalization found with the name:", normalization))
    }
  }
  else
  {
    normVals <- rep(1, nrow(peakMatrix$pos))
  }
  
  offsetPosMat <- OffsetRasterPosMats(peakMatrix$pos, peakMatrix$numPixels, offset_horizontal)
  rMSI::PlotValues(offsetPosMat, peakMatrix[[matrix]][,column]/normVals, rotate = rotate, scale_title = sprintf("m/z %.4f", peakMatrix$mass[column]), pixel_size_um = pixel_size_um)
}

#' plotClusterImage.
#'
#' Plot a segmentation image with the user-given clusters.
#'
#' @param peakMatrix the peak matrix in an rMSIproc object.
#' @param clusters a vector with integer number according the cluster of each pixel.
#' @param offset_horizontal if true images will be laydout horizontally.
#' @param rotate rotation to apply.
#' @param pixel_size_um the pixel resolution in um.
#'
#' @return a vector with the used color for each cluster sorted according clustering numering in assending order.
#'
plotClusterImage <- function(peakMatrix, clusters,  offset_horizontal = T, rotate = 0, pixel_size_um = 100)
{
  offsetPosMat <- OffsetRasterPosMats(peakMatrix$pos, peakMatrix$numPixels, offset_horizontal)
  return(rMSI::PlotClusterImage(offsetPosMat, clusters, rotate, pixel_size_um))
}
