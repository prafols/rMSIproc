#########################################################################
#     
#     Copyright (C) 2018 Lluc Sementé Fernàndez
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
#' Deisotoping.
#' 
#' Finds and evaluate isotope candidates for each ion mass.
#' 
#'
#' @param PeakMtx An rMSIprocPeakMatrix.
#' @param IsoNumber Number of isotopes to be found.
#' @param Tolerance Mass tolerance for the candidates in scans or ppms.
#' @param ScoreThreshold Score value to consider a ion mass a good isotope candidate. Only the ions that have this number or greater will undergo the following isotope searching stages.
#' @param ImageVector The mass channels vector of the imaging dataset.
#' 
#' @return A list containing the results of the scores and other kind of information.
#' 
#' @export

Deisotoping <- function(PeakMtx, IsoNumber = 2,Tolerance = 10, ScoreThreshold = 0.8,ToleranceUnits = "ppm", ImageVector = NULL)
{
  if(ToleranceUnits == "ppm")
  {
    InScans <- FALSE
  } else if(ToleranceUnits == "scan")
    {
      InScans <- TRUE
    } else
      {
        writeLines("Error: ToleranceUnits must be 'ppm' or 'scan'")
        break
      }
  
  if((ToleranceUnits == "scan") & (length(PeakMtx$mass) > length(ImageVector)))
  {
    writeLines("Error: If ToleranceUnits is 'scan' then ImageVector must be the mass channels vector of the rMSI image. (rMSIObj$mass)")
    break
  }

  ############################## Calling the C++ method ##############################
  r <- list()  
  r <- IsotopeAnnotator(length(PeakMtx$mass),     #number of mass peaks
                        length(ImageVector),      #number of massChannels
                        sum(PeakMtx$numPixels),   #number of pixels
                        IsoNumber,                #number of isotopes
                        PeakMtx$intensity,        #peak matrix
                        PeakMtx$mass,             #matrix mass axis
                        ImageVector,              #image mass axis 
                        Tolerance,                #tolerance in ppm or scans
                        ScoreThreshold,           #score threshold
                        InScans)                  #tolerance units
  
  ############################## Giving format to the output ##############################
  NullVec <- matrix(nrow = length(r)-1, ncol = length(r[[1]]))
  Nullcnt <- 0

  #Naming the list
  for(i in 2:length(r))
  {
    names(r[[i]]) <- signif(r[[1]],digits = 8)
    Nullcnt <- 0
    for(j in 1:length(r[[i]]))
    {
      if(!is.null(r[[i]][[j]]))
      {
        row.names(r[[i]][[j]]) <- c("M+N mass", "M+0 index", "M+N index","ppm error", "Final score","Morphology score", "Intensity score", "Model slope", "Number of C atoms")
        colnames(r[[i]][[j]]) <- (1:(dim(r[[i]][[j]])[2]))
      }else
      {
        Nullcnt <- Nullcnt + 1
        NullVec[i-1,Nullcnt] <- j
      }
    }
    r[[i]] <- r[[i]][-(NullVec[i-1,1:Nullcnt])]
  }
 
  #Removing the mass vector and the empty isotope lists
  r <- r[-1]
  flagRemList <- FALSE
  RemList <- c()
  
  for(i in 1:length(r))
  {
    if(length(r[[i]])==0)
    {
      flagRemList <- TRUE
      RemList <- c(RemList,i)
    }
  }
  
  if(flagRemList)
  {
    r <- r[-RemList]
  }
  
  #Creating the isotope index vector
  CompVec <- c()
  for(i in (1:(length(r))))
  {
    for(j in (1:length(r[[i]])))
    {
      if(!is.null(r[[i]][[j]]))
      {
        if(r[[i]][[j]][5,which.max(r[[i]][[j]][5,])] >= ScoreThreshold)
        {
          CompVec <- c(CompVec, r[[1]][[j]][3,which.min(abs(r[[1]][[j]][4,which(r[[1]][[j]][5,]>= ScoreThreshold)]))])  
        }
      }  
    }
  }
  names(r) <- paste("M+",1:length(r),sep = "")
  r$isotopeIndex <- sort(unique(CompVec))
  
  
  MonoVec <- c()
  for(i in 1:length(r[[1]]))
  {
    if(!is.null(r[[1]][[i]]))
    {
      if(r[[1]][[i]][5,which.min(abs(r[[1]][[i]][4,]))] >= ScoreThreshold)
      {
        MonoVec <- c(MonoVec, r[[1]][[i]][2,which.min(abs(r[[1]][[i]][4,]))])
      }
    }
  }
  r$MonoisotopeIndex <- sort(unique(MonoVec))
  
  return(r)
}