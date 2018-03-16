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
#include "mtaverage.h"
using namespace Rcpp;

MTAverage::MTAverage(AvrgDef imgRunInfo) : 
  ThreadingMsiProc(imgRunInfo.numOfThreads, false,
                   imgRunInfo.fileNames, imgRunInfo.massChannels,
                   imgRunInfo.numRows, imgRunInfo.dataType)
{
  
  massCh = imgRunInfo.massChannels;
  numCubes = imgRunInfo.fileNames.length();
  numPixels = 0;
  for (int i = 0; i < numCubes; i++)
  {
  numPixels += imgRunInfo.numRows[i];
  }
  sm = new double*[numCubes];
  for(int i = 0; i < numCubes ; i++)
  {
    sm[i] = new double[imgRunInfo.massChannels];
  }
  
  for (int i = 0; i < numCubes ; i++)
  {
    for (int j = 0; j < massCh ; j++)
    {
    sm[i][j] = 0;
    }
  }
}

MTAverage::~MTAverage()
{
  for (int i = 0; i < numCubes; i++)
  {
     delete[] sm[i];
  }
  
  delete[] sm;
}

NumericVector MTAverage::Run()
{
  NumericVector avg(massCh);
  for (int j = 0; j < massCh; j++)
  {
    avg[j] = 0;
  }
  //Run smoothing in multi-threading
  runMSIProcessingCpp();
  
  //Merging all the cube's partial average spectrum and getting the total average spectrum
  for (int i = 0; i < numCubes; i++)
  {
    for (int j = 0; j < massCh; j++)
    {
      avg[j] += sm[i][j];
    }
  }
  
  for (int j = 0; j < massCh; j++)
  {
    avg[j] /= numPixels;    //Dividir per el nombre total de pixels.
  }
  
  return avg;
}



void MTAverage::ProcessingFunction(int threadSlot)
{
  //Perform the average value of each mass channel in the current loaded cube
  for (int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
       sm[cubes[threadSlot]->cubeID][k] += cubes[threadSlot]->data[j][k];
    }
  }
}


// [[Rcpp::export]]
NumericVector AverageSpectrumC(StringVector fileNames, 
                          int massChannels, IntegerVector numRows,
                          String dataType, int numOfThreads)
{
  //Copy R data to C arrays
  int numRowsC[fileNames.length()];
  memcpy(numRowsC, numRows.begin(), sizeof(int)*fileNames.length());
  
  MTAverage::AvrgDef myProcParams;
  myProcParams.dataType = dataType;
  myProcParams.fileNames = fileNames;
  myProcParams.massChannels = massChannels;
  myProcParams.numOfThreads = numOfThreads;
  myProcParams.numRows = numRowsC; 
  
  MTAverage myAverage(myProcParams);
  return myAverage.Run();
}
