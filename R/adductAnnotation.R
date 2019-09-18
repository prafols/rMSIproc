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
  adducts <- list()
  
  #namber of adduct combinations
  combinations <- 0
  for(i in 1:nrow(adductDataFrame)-1)
  {
    combinations <- combinations + i
  }
  
  #adduct labels for the output
  name1 = 1
  name2 = 2
  lim = nrow(adductDataFrame)
  namesVector <- c()
  firstnameVector <- c()
  secondnameVector <- c()
  for(i in 1:combinations)
  {
    namesVector[i] <- paste(adductDataFrame$name[name1]," & ",adductDataFrame$name[name2],sep = "")
    firstnameVector[i] <- as.character(adductDataFrame$name[name1])
    secondnameVector[i] <- as.character(adductDataFrame$name[name2])
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

  
  ##### Calling the C++ method #####
  adducts <- C_adductAnnotation(length(isotopeObj$monoisotopicPeaks),
                                nrow(adductDataFrame), 
                                tolerance,
                                length(PeakMtx$mass),
                                PeakMtx$mass[sort(isotopeObj$monoisotopicPeaks)],
                                adductDataFrame$mass,
                                isotopeObj$M1,
                                ord,
                                PeakMtx$mass,
                                PeakMtx$intensity,
                                sum(PeakMtx$numPixels))   
  #return(adducts)
  ##### Output format #####
  if (length(adducts) >= 2) 
  {
    tmp <- adducts[[length(adducts)]]
    adducts <- adducts[-length(adducts)]
    names(adducts) <- format(tmp,nsmall = 4)
  
    for(i in 1:length(adducts))
    {
      name1 <- firstnameVector[(adducts[[i]][1]+1)]
      name2 <- secondnameVector[(adducts[[i]][1]+1)]
      mainname <- namesVector[adducts[[i]][1]+1]
      adducts[[i]] <- as.matrix(adducts[[i]][-1,])
  
      row.names(adducts[[i]]) <- c(paste("M+",name1," mass",sep = ""),
                                   paste("M+",name1," index",sep = ""),
                                   paste("M+",name2," mass",sep = ""),
                                   paste("M+",name2," index",sep = ""),
                                   "C atoms mean",
                                   "C atoms avg. error",
                                   "Correlation",
                                   "ppm error")
       colnames(adducts[[i]]) <- mainname
    }
    adducts <- adducts[order(as.numeric(names(adducts)))]
    return(adducts)
  }else
    {
      return(NULL)
    }
}
