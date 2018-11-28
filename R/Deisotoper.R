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

Deisotoping <- function(PeakMtx, IsoNumber, Scans, Score, MassChannelsVec)
{
  
  #Calling the C++ method
  r <- list()  
  r <- IsotopeAnnotator(length(PeakMtx$mass), length(MassChannelsVec), 
                     sum(PeakMtx$numPixels), IsoNumber, PeakMtx$intensity, 
                     PeakMtx$mass, MassChannelsVec, Scans, Score)
  
  #Giving format to the output
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
      row.names(r[[i]][[j]]) <- c("Candidate mass", "Final score", "Morphology score", "Intensity score", "Model slope", "Number of C atoms")
      colnames(r[[i]][[j]]) <- (1:(dim(r[[i]][[j]])[2]))
      }else
       {
        Nullcnt <- Nullcnt + 1
        NullVec[i-1,Nullcnt] <- j
       }
    }
    r[[i]] <- r[[i]][-(NullVec[i-1,1:Nullcnt])]
  }

  #Removing the mass vector
  r <- r[-1]

  return(r)
}