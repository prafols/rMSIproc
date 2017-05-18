/*************************************************************************
 *     rMSIproc - R package for MSI data processing
 *     Copyright (C) 2017 Pere Rafols Soler
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/

#ifndef MT_REPLACEZEROPEAKS_H
#define MT_REPLACEZEROPEAKS_H
#include <Rcpp.h>
#include "peakpicking.h"
#include "threadingmsiproc.h"

class MTReplaceZeroPeaks : public ThreadingMsiProc 
{
  public:
    //Data structure used to completely define the processing pipeline
    typedef struct
    {
      int peakWinSize; //Windows size used for peak-picking
      double *massAxis; //Array containing the mass axis
      int massChannels; //Number of points in the mass axis array
      int peakInterpolationUpSampling; //Upsampling value for peak interpolation
  
      Rcpp::String basePath; //Full path where ramdisks are stored
      Rcpp::StringVector fileNames; //Filname of each ramdisk file
      int *numRows; //An array containing the number of rows stored in each ramdisk file. The length ot this array is the length of fileNames
      Rcpp::String dataType; //A string with the data type
  
      int numOfThreads;
    }ImgProcDef;  
    
    MTReplaceZeroPeaks(ImgProcDef imgRunInfo, Rcpp::List peakMatrix);
    ~MTReplaceZeroPeaks();
    
    //Exectur a full imatge processing using threaded methods
    Rcpp::List Run(); 
    
  private:
    int *replacedZerosCounters;
    int *mass_index;
    Rcpp::NumericVector *pkMatmass;
    Rcpp::NumericMatrix *pkMatintensity;
    Rcpp::NumericMatrix *pkMatarea;
    Rcpp::NumericMatrix *pkMatsnr;
    PeakPicking **peakObj;
  
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
};
#endif
