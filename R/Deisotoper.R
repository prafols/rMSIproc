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
#' @param IsoNumber Integer. Number of isotopes to be found.
#' @param Tolerance Integer. Mass tolerance for the candidates in scans or ppms.
#' @param ScoreThreshold Numeric. Score value to consider a ion mass a good isotope candidate. Only the ions that have this number or greater will undergo the following isotope searching stages.
#' @param ToleranceUnits String. Must be "ppm" or "scan". If ToleranceUnits is 'scan' then ImageVector must be the mass channels vector of the rMSI image (rMSIObj$mass), else tolerance will be interpreted in ppms 
#' @param ImageVector Numeric Vector. The mass channels vector of the imaging dataset.
#' 
#' @return A list containing the results of the scores and other kind of information.
#' 
#' @export

Deisotoping <- function(PeakMtx, IsoNumber = 2, Tolerance = 10, ScoreThreshold = 0.8, ToleranceUnits = "ppm", ImageVector = NULL)
{
  ##### Type check #####
  if(ToleranceUnits == "ppm")
  {
    InScans <- FALSE
  } else if(ToleranceUnits == "scan")
    {
      InScans <- TRUE
    } else
      {
        stop("ToleranceUnits must be 'ppm' or 'scan'")
      }
  
  if((ToleranceUnits == "scan") & (length(PeakMtx$mass) > length(ImageVector)))
  {
    stop("If ToleranceUnits is 'scan' then ImageVector must be the mass channels vector of the rMSI image. (rMSIObj$mass)")
  }

  if(!(IsoNumber%%1 == 0))
  {
    stop("Error: IsoNumber must be an integer")
  }
  
  if((ToleranceUnits == "scan") & !(Tolerance%%1 == 0))
  {
    stop("Error: If ToleranceUnits is 'scan' then Tolerance must be integer")
  }
  
  ##### Calling the C++ method #####
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
  

  ##### Names the results list #####
  NullVec <- matrix(nrow = length(r)-1, ncol = length(r[[1]]))
  Nullcnt <- 0

  for(i in 2:length(r))
  {
    names(r[[i]]) <- signif(r[[1]],digits = 8)
    Nullcnt <- 0
    for(j in 1:length(r[[i]]))
    {
      if(!is.null(r[[i]][[j]]))
      {
        row.names(r[[i]][[j]]) <- c("M+N mass","ppm error","M+0 index", "M+N index","Final score","Morphology score","Intensity score", "Mass error score", "Isotopic ratio", "Number of C atoms")
        colnames(r[[i]][[j]]) <- (1:(dim(r[[i]][[j]])[2]))
      }else
        {
          Nullcnt <- Nullcnt + 1
          NullVec[i-1,Nullcnt] <- j
        }
    }
    r[[i]] <- r[[i]][-(NullVec[i-1,1:Nullcnt])]
  }
 
  
  ##### Removes the mass vector and the empty isotope lists #####
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
  

  ##### Creats the isotopic & monoisotopic ions index vectors #####
  MonoVec <- c()
  CompVec <- c()
  for(i in (1:(length(r))))
  {
    for(j in (1:length(r[[i]])))
    {
      if((!is.null(r[[i]][[j]])) & (r[[i]][[j]][5,which.max(r[[i]][[j]][5,])] >= ScoreThreshold))
      {
        if(i == 1)
        {
          MonoVec <- c(MonoVec, r[[1]][[j]][3,which.max(r[[1]][[j]][5,])])
        }
        CompVec <- c(CompVec, r[[i]][[j]][4,which.max(r[[i]][[j]][5,])])
      }
    }
  }
  names(r) <- paste("M",1:length(r),sep = "")
  r$isotopicPeaks <- sort(unique(CompVec))
  r$monoisotopicPeaks <- sort(unique(MonoVec))
  
  
  return(r)
}