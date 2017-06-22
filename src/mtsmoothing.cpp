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

#include <Rcpp.h>
#include <cmath>
#include "mtsmoothing.h"
using namespace Rcpp;

MTSmoothing::MTSmoothing(ImgProcDef imgRunInfo) : 
  ThreadingMsiProc( imgRunInfo.numOfThreads, true, imgRunInfo.fileNames, imgRunInfo.massChannels, imgRunInfo.numRows, imgRunInfo.dataType )
{
  smObj = new Smoothing*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    smObj[i] = new Smoothing(imgRunInfo.SmoothingKernelSize);
  }
}

MTSmoothing::~MTSmoothing()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete smObj[i];
  }
  delete[] smObj;
}

void MTSmoothing::Run()
{
  //Run smoothing in mutli-threading
  runMSIProcessingCpp();
}

void MTSmoothing::ProcessingFunction(int threadSlot)
{
  //Perform Savitzky-Golay smoothing of each spectrum in the current loaded cube
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    smObj[threadSlot]->smoothSavitzkyGolay(cubes[threadSlot]->data[j], cubes[threadSlot]->ncols);
  }
}

// [[Rcpp::export]]
void FullImageSmoothing ( StringVector fileNames, 
                                int massChannels, IntegerVector numRows,
                                String dataType, int numOfThreads, 
                                int SmoothingKernelSize = 5 )
{
  //Copy R data to C arrays
  int numRowsC[fileNames.length()];
  memcpy(numRowsC, numRows.begin(), sizeof(int)*fileNames.length());
  
  MTSmoothing::ImgProcDef myProcParams;
  myProcParams.dataType = dataType;
  myProcParams.fileNames = fileNames;
  myProcParams.massChannels = massChannels;
  myProcParams.numOfThreads = numOfThreads;
  myProcParams.numRows = numRowsC; 
  myProcParams.SmoothingKernelSize = SmoothingKernelSize;
  
  MTSmoothing mySmoothing(myProcParams);
  mySmoothing.Run();
}
