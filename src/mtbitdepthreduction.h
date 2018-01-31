/*************************************************************************
 *     rMSIproc - R package for MSI data processing
 *     Copyright (C) 2014 Pere Rafols Soler
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

#ifndef MT_BITDEPTHREDUCTION_H
  #define MT_BITDEPTHREDUCTION_H
#include <Rcpp.h>
#include "threadingmsiproc.h"
#include "noiseestimation.h"

class MTBitDepthReduction : public ThreadingMsiProc 
{
  public:
    //Data structur used to completly define the processing pipeline
    typedef struct
    {
      int massChannels; //Number of points in the mass axis array
      Rcpp::StringVector fileNames; //Filname of each ramdisk file
      int *numRows; //An array containing the number of rows stored in each ramdisk file. The length ot this array is the length of fileNames
      Rcpp::String dataType; //A string with the data type
      int numOfThreads;
      int noiseWinSize; // The windows size used to estimate noise
    }ImgProcDef;  
    
    MTBitDepthReduction(ImgProcDef imgRunInfo);
    ~MTBitDepthReduction();
    
    //Exectur a full imatge processing using threaded methods
    void Run();
    
    //Single spectrum bit depth reduction
    void BitDepthReduction(double *data, int dataLength, int noiseModelThreadSlot);
    
  private:
    NoiseEstimation **noiseModel;
    int NoiseWinSize;
    unsigned long long maskLUT_double[52];
    
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
};
#endif
