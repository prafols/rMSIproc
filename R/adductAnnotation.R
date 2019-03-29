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
#' Finds and evaluate isotope candidates for each ion mass.
#' 
#' @param monoisitopeMassVector Integer. Number of isotopes to be found.
#' @param adductDataFrame Integer. Mass tolerance for the candidates in scans or ppms.
#' @param scoreThreshold Numeric. Score value to consider a ion mass a good isotope candidate. Only the ions that have this number or greater will undergo the following isotope searching stages.
#
#' @return A list containing the results of the test and other kind of information.
#' 

adductAnnotation <- function(monoisitopeMassVector, adductDataFrame, tolerance)
{
  result <- list()  
  
  ##### Calling the C++ method #####
  result <- adductAnnotation(length(monoisitopeMassVector),
                             nrow(adductDataFrame), 
                             tolerance, 
                             monoisitopeMassVector,
                             adductDataFrame$mass)   
 
  ##### Output format #####
  
  return(result)
}