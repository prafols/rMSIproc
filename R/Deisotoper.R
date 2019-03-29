#########################################################################
#     
#     Copyright (C) 2018 Lluc Semente Fernandez
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
#' Deisotoping
#' 
#' Finds and evaluate isotope candidates for each ion mass.
#' 
#'
#' @param PeakMtx List. An rMSIprocPeakMatrix. Must contain at least the following categories: \itemize{
#' \item PeakMtx$intensity. A matrix containg the intensities of the peaks for each pixel (rows = pixels, cols = peaks).  
#' \item PeakMtx$mass. A vector containg the masses of each peak. Must be in the same order with the columns of the intensity marix.   
#' \item PeakMtx$numPixels. Number of pixels (rows in your matrix).   
#' }
#' @param isoNumber Integer. Number of isotopes to be found.
#' @param tolerance Integer. Mass tolerance for the candidates in scans or ppms.
#' @param scoreThreshold Numeric. Score value to consider a ion mass a good isotope candidate. Only the ions that have this number or greater will undergo the following isotope searching stages.
#' @param toleranceUnits String. Must be 'ppm' or 'scan'. If ToleranceUnits is 'scan' then ImageVector must be the mass channels vector of the rMSI image (rMSIObj$mass).
#' @param imageVector Numeric Vector. The mass channels vector of the imaging dataset containing all the scans.
#' 
#' @return A list containing the results of the test and other kind of information.
#' 

Deisotoping <- function(PeakMtx, isoNumber = 2, tolerance = 30, scoreThreshold = 0.8, toleranceUnits = "ppm", imageVector = NULL)
{
  result <- list()  
  
  ##### Type check #####
  {
    if(!(is.matrix(PeakMtx$intensity) & (length(PeakMtx$mass)>0)))
    {
      stop("Error: PeakMtx has not the rMSIprocPeakMatrix format. Look at help for more information.")
    }
    
    
    if(isoNumber < 1)
    {
      stop("Error: IsoNumber must be 1 or a greater integer")
    }
    
    if ((scoreThreshold > 1) | (scoreThreshold < 0))
    {
      stop("Error: ScoreThreshold must be between 0 and 1")
    }
    
    if(toleranceUnits == "ppm")
    {
      InScans <- FALSE
      imageVector <- c(1,2)
    } else if(toleranceUnits == "scan")
    {
      InScans <- TRUE
    } else
    {
      stop("Error: ToleranceUnits must be 'ppm' or 'scan'")
    }
    
    if((toleranceUnits == "scan") & (length(PeakMtx$mass) > length(imageVector)))
    {
      stop("Error: If ToleranceUnits is 'scan' then ImageVector must be the mass channels vector of the rMSI image. (rMSIObj$mass)")
    }
    
    if((toleranceUnits == "ppm") & (tolerance>500))
    {
      writeLines("Maximum allowed tolerance is 500. Replacing the value")
      tolerance = 500
    }
    
    if(!(isoNumber%%1 == 0))
    {
      stop("Error: IsoNumber must be an integer")
    }
    
    if((toleranceUnits == "scan") & !(tolerance%%1 == 0))
    {
      stop("Error: If ToleranceUnits is 'scan' then Tolerance must be integer")
    }
    
    if((toleranceUnits == "scan") & (tolerance>length(imageVector)))
    {
      stop("Error: Not enough scans in your mass axis, reduce the tolerance")
    }
  }
  ##### Calling the C++ method #####
  result <- isotopeAnnotator(length(PeakMtx$mass),   #number of mass peaks
                        length(imageVector),         #number of massChannels
                        sum(PeakMtx$numPixels),      #number of pixels
                        isoNumber,                   #number of isotopes
                        PeakMtx$intensity,           #peak matrix
                        PeakMtx$mass,                #matrix mass axis
                        imageVector,                 #image mass axis 
                        tolerance,                   #tolerance in ppm or scans
                        scoreThreshold,              #score threshold
                        InScans)                     #tolerance units
 
  ##### Output format #####
  result <- DeisotopingOutputFormat(result,scoreThreshold)
  
  return(result)
}

#' DeisotopingOutputFormat
#' 
#' Gives format to the output. 

DeisotopingOutputFormat <- function(r,ScoreThreshold)
{
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
        row.names(r[[i]][[j]]) <- c("M+N mass","Abs ppm error","M+0 index", "M+N index","Final score","Morphology score","Intensity score", "Mass error score", "Isotopic ratio", "Number of C atoms")
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
  writeLines(paste("Number of monoisotopic ions found =",length(r$monoisotopicPeaks)))
  
  return(r)
}