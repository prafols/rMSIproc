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
#' 
#' 
#' @param isotopeObj List. Results from the isotopes test.
#' @param massVector  Vector. Mass vector from wich ceck adduct relationship.
#' @param adductDataFrame Data frame with two columns. $name, with the names of adducts to be found, $mass, with each masses.
#' @param tolerance Integer. Mass error tolerance.
#' @return 
#' 

adductAnnotation <- function(isotopeObj, massVector, adductDataFrame, tolerance)
{
  adducts <- list()
  result <- list()
  result$isotopes <- isotopeObj 
  
  combinations <- 0
  for(i in 1:nrow(adductDataFrame))
  {
    combinations <- combinations + 1
  }
  
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
  
  ##### Calling the C++ method #####
  adducts <- C_adductAnnotation(length(isotopeObj$monoisotopicPeaks),
                                            nrow(adductDataFrame), 
                                            tolerance,
                                            length(massVector),
                                            sort(massVector[isotopeObj$monoisotopicPeaks])-1,
                                            adductDataFrame$mass,
                                            isotopeObj$M1,
                                            order(isotopeObj$monoisotopicPeaks)-1,
                                            massVector)   
  ##### Output format #####

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
                                 paste("M+",name1," slope",sep = ""),
                                 paste("M+",name2," mass",sep = ""),
                                 paste("M+",name2," index",sep = ""),
                                 paste("M+",name2," slope",sep = ""),
                                 "Correlation",
                                 "ppm error")
     colnames(adducts[[i]]) <- mainname
  }

  result$adducts <- adducts
  
  
  return(result)
}