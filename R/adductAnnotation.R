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
#' adductAnnotation
#' 
#' Given the monoisotopic ions found by the isotopes test, founds the possible adduct pairs considering the elements in the adductDataFrame.
#' 
#' @param isotopeObj List. Results from the isotopes test.
#' @param PeakMtx List. An rMSIprocPeakMatrix. Must contain at least the following categories: \itemize{
#' \item PeakMtx$intensity. A matrix containg the intensities of the peaks for each pixel (rows = pixels, cols = peaks).  
#' \item PeakMtx$mass. A vector containg the masses of each peak. Must be in the same order with the columns of the intensity marix.   
#' \item PeakMtx$numPixels. Number of pixels (rows in your matrix).   
#' }
#' @param adductDataFrame Data frame with two columns. $name, with the names of adducts to be found, $mass, with each masses.
#' @param tolerance Integer. Mass error tolerance in ppm.
#' 
#' @return A list containing the results of the test and other kind of information.
#' 

adductAnnotation <- function(isotopeObj, PeakMtx, adductDataFrame, tolerance)
{
  if(!is.data.frame(adductDataFrame))
  { 
    stop("add.adducts must be a data frame")
  }
  
  if(tolerance>500)
  {
    writeLines("Maximum allowed tolerance is 500. Replacing the value.")
    tolerance <- 500
  }
  
  adducts <- list()
  adductDataFrame <- adductDataFrame[order(adductDataFrame$mass,decreasing = T),]
  #namber of adduct combinations
  combinations <- 0
  for(i in 1:nrow(adductDataFrame)-1)
  {
    combinations <- combinations + i
  }
  
  #adduct labels for the output
  name1 <- 1
  name2 <- 2
  lim <- nrow(adductDataFrame)
  namesVector <- c()
  firstnameVector <- c()
  secondnameVector <- c()
  firstnamePriority <- c()
  secondnamePriority <- c()
  for(i in 1:combinations)
  {
    namesVector[i] <- paste(adductDataFrame$name[name1]," & ",adductDataFrame$name[name2],sep = "")
    firstnameVector[i] <- as.character(adductDataFrame$name[name1])
    secondnameVector[i] <- as.character(adductDataFrame$name[name2])
    firstnamePriority[i] <- adductDataFrame$priority[name1]
    secondnamePriority[i] <- adductDataFrame$priority[name2]
    name2 <- name2 + 1
    if(name2 > lim)
    {
      name1 <- name1 + 1
      name2 <- name1 + 1
    }
  }
  
  #monoisotopic to list order
  ord <- c()
  for (i in 1:length(isotopeObj$monoisotopicPeaks)) 
  {
    ord <- c(ord,which.min(abs(PeakMtx$mass[sort(isotopeObj$monoisotopicPeaks)[i]]-as.numeric(names(isotopeObj$M1)))))
  }
  ord <- ord-1
  
  #label axis
  labelAxis <- rep(0, times = length(PeakMtx$mass))
  labelAxis[isotopeObj$monoisotopicPeaks] <- 1
  labelAxis[isotopeObj$isotopicPeaks] <- 2
  
  M1isotopes <- isotopeObj$M1
  
  for(i in 1:length(M1isotopes))
  {
    M1isotopes[[i]] <- M1isotopes[[i]][,which.max(M1isotopes[[i]][5,])]
  }
  
  ##### Calling the C++ method #####
  adducts <- C_adductAnnotation(length(isotopeObj$monoisotopicPeaks),
                                nrow(adductDataFrame), 
                                tolerance,
                                length(PeakMtx$mass),
                                PeakMtx$mass[sort(isotopeObj$monoisotopicPeaks)],
                                adductDataFrame$mass,
                                M1isotopes,
                                ord,
                                PeakMtx$mass,
                                PeakMtx$intensity,
                                sum(PeakMtx$numPixels),
                                labelAxis,
                                sort(isotopeObj$monoisotopicPeaks)-1)   

  for(i in 1:2)
  {
    adducts[[i]] <- adducts[[i]][-((length(adducts[[i]])-sum(unlist(lapply(adducts[[i]], is.null)))+1):length(adducts[[i]]))]
  }
  names(adducts) <- c("A","B")
  
  ## Data frame for A quality adducts ##
  if(length(adducts$A) > 0)
  {
    adductsA <- data.frame(
    NeutralMass = rep(0, times = length(adducts$A)),
    Adducts = rep(0, times = length(adducts$A)),
    Adduct1Mass = rep(0, times = length(adducts$A)),
    Adduct2Mass = rep(0, times = length(adducts$A)),
    IsotopeIntensityRatioMean = rep(0, times = length(adducts$A)),
    IsotopeIntensityRatioStdError = rep(0, times = length(adducts$A)),
    Correlation = rep(0, times = length(adducts$A)),
    MassError = rep(0, times = length(adducts$A)),
    Adduct1Index = rep(0, times = length(adducts$A)),
    Adduct2Index = rep(0, times = length(adducts$A)),
    AdductPriority = rep(0, times = length(adducts$A)))
    for(i in 1:length(adducts$A))
    {
      name1 <- firstnameVector[(adducts$A[[i]][1]+1)]
      name2 <- secondnameVector[(adducts$A[[i]][1]+1)]
      if(adducts$A[[i]][3] > adducts$A[[i]][5])
      {
        adductsA$Adducts[i] <- paste("[M",name1,"] & [M",name2,"]",sep="")
      }
      else
      {
        adductsA$Adducts[i] <- paste("[M",name2,"] & [M",name1,"]",sep="")
      }
      adductsA$NeutralMass[i]                   <- adducts$A[[i]][2]
      adductsA$Adduct1Mass[i]                   <- adducts$A[[i]][3]
      adductsA$Adduct1Index[i]                  <- adducts$A[[i]][4]
      adductsA$Adduct2Mass[i]                   <- adducts$A[[i]][5]
      adductsA$Adduct2Index[i]                  <- adducts$A[[i]][6]
      adductsA$IsotopeIntensityRatioMean[i]     <- adducts$A[[i]][7]
      adductsA$IsotopeIntensityRatioStdError[i] <- adducts$A[[i]][8]
      adductsA$Correlation[i]                   <- adducts$A[[i]][9]
      adductsA$MassError[i]                     <- adducts$A[[i]][10] 
      adductsA$AdductPriority[i]                <- (firstnamePriority[(adducts$A[[i]][1]+1)]) + (secondnamePriority[(adducts$A[[i]][1]+1)])
    }
    adductsA$Correlation <- signif(adductsA$Correlation, digits = 3)
    adductsA$MassError <- trunc(adductsA$MassError) + signif(adductsA$MassError - trunc(adductsA$MassError), digits = 3)
    adductsA$Adduct1Mass <- trunc(adductsA$Adduct1Mass) + signif(adductsA$Adduct1Mass - trunc(adductsA$Adduct1Mass), digits = 4)
    adductsA$Adduct2Mass <- trunc(adductsA$Adduct2Mass) + signif(adductsA$Adduct2Mass - trunc(adductsA$Adduct2Mass), digits = 4)
    adductsA$NeutralMass <- trunc(adductsA$NeutralMass) + signif(adductsA$NeutralMass - trunc(adductsA$NeutralMass), digits = 4)
    adductsA$IsotopeIntensityRatioMean <- trunc(adductsA$IsotopeIntensityRatioMean) + signif(adductsA$IsotopeIntensityRatioMean - trunc(adductsA$IsotopeIntensityRatioMean), digits = 3)
    adductsA$IsotopeIntensityRatioStdError <- trunc(adductsA$IsotopeIntensityRatioStdError) + signif(adductsA$IsotopeIntensityRatioStdError - trunc(adductsA$IsotopeIntensityRatioStdError), digits = 4)
    
    adductsA <- adductsA[order(adductsA$NeutralMass,decreasing = F),]
    row.names(adductsA) <- NULL
    adducts$A <- adductsA[order(adductsA$NeutralMass,decreasing = F),]
  }
  else {adducts <- adducts[-1]}
    
  ## Data frame for B quality adducts ##
  if(length(adducts$B) > 0 )
  {
    adductsB <- data.frame(
    NeutralMass = rep(0, times = length(adducts$B)),
    Adducts = rep(0, times = length(adducts$B)),
    Adduct1Mass = rep(0, times = length(adducts$B)),
    Adduct2Mass = rep(0, times = length(adducts$B)),
    Correlation = rep(0, times = length(adducts$B)),
    MassError = rep(0, times = length(adducts$B)),
    Adduct1Index = rep(0, times = length(adducts$B)),
    Adduct2Index = rep(0, times = length(adducts$B)))
    for(i in 1:length(adducts$B))
    {
      name1 <- firstnameVector[(adducts$B[[i]][1]+1)]
      name2 <- secondnameVector[(adducts$B[[i]][1]+1)]
      if(adducts$B[[i]][3] > adducts$B[[i]][5])
      {
        adductsB$Adducts[i] <- paste("[M",name1,"] & [M",name2,"]",sep="")
      }
      else
      {
        adductsB$Adducts[i] <- paste("[M",name2,"] & [M",name1,"]",sep="")
      }
      adductsB$NeutralMass[i]    <- adducts$B[[i]][2]
      adductsB$Adduct1Mass[i]    <- adducts$B[[i]][3]
      adductsB$Adduct1Index[i]   <- adducts$B[[i]][4]
      adductsB$Adduct2Mass[i]    <- adducts$B[[i]][5]
      adductsB$Adduct2Index[i]   <- adducts$B[[i]][6]
      adductsB$Correlation[i]    <- adducts$B[[i]][7]
      adductsB$MassError[i]      <- adducts$B[[i]][8]
      adductsB$AdductPriority[i] <- firstnamePriority[(adducts$B[[i]][1]+1)] + secondnamePriority[(adducts$B[[i]][1]+1)]
    }
    adductsB$Correlation <- signif(adductsB$Correlation, digits = 3)
    adductsB$MassError <- trunc(adductsB$MassError) + signif(adductsB$MassError - trunc(adductsB$MassError), digits = 3)
    adductsB$Adduct1Mass <- trunc(adductsB$Adduct1Mass) + signif(adductsB$Adduct1Mass - trunc(adductsB$Adduct1Mass), digits = 4)
    adductsB$Adduct2Mass <- trunc(adductsB$Adduct2Mass) + signif(adductsB$Adduct2Mass - trunc(adductsB$Adduct2Mass), digits = 4)
    adductsB$NeutralMass <- trunc(adductsB$NeutralMass) + signif(adductsB$NeutralMass - trunc(adductsB$NeutralMass), digits = 4)
    
    adductsB <- adductsB[order(adductsB$NeutralMass,decreasing = F),]
    row.names(adductsB) <- NULL
    adducts$B <- adductsB[order(adductsB$NeutralMass,decreasing = F),]
  }
  else {adducts <- adducts[-2]}

  return(adducts)
}
