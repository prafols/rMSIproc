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
#include "mtbitdepthreduction.h"
using namespace Rcpp;

#define MIN_MANTISSA_BITS 4
#define NOISE_THRESHOLD_LOWER 0.5
#define NOISE_THRESHOLD_UPPER 10.0

MTBitDepthReduction::MTBitDepthReduction(ImgProcDef imgRunInfo) : 
  ThreadingMsiProc( imgRunInfo.numOfThreads, true, imgRunInfo.fileNames, imgRunInfo.massChannels, imgRunInfo.numRows, imgRunInfo.dataType )
{
  NoiseWinSize = imgRunInfo.noiseWinSize;
  noiseModel = new NoiseEstimation*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    noiseModel[i] = new NoiseEstimation(imgRunInfo.massChannels);
  }
  
  //Fill the LUT with all possible resolutions from 1 to 52 bits
  maskLUT_double[0] = 0xFFF8000000000000;
  for( int i = 1; i < 52; i++)
  {
    maskLUT_double[i] = (maskLUT_double[i-1] >> 1) | 0x8000000000000000; 
  }
}

MTBitDepthReduction::~MTBitDepthReduction()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete noiseModel[i];
  }
  delete[] noiseModel;
}

void MTBitDepthReduction::Run()
{
  //Run bit depth reduction in mutli-threading
  runMSIProcessingCpp();
}

void MTBitDepthReduction::ProcessingFunction(int threadSlot)
{
  //Perform Savitzky-Golay smoothing of each spectrum in the current loaded cube
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    BitDepthReduction(cubes[threadSlot]->data[j], cubes[threadSlot]->ncols, threadSlot);
  }
}

void MTBitDepthReduction::BitDepthReduction(double *data, int dataLength, int noiseModelThreadSlot)
{
  
  int resolution_bits;
  double *noise_floor = new double[dataLength];
  memcpy(noise_floor, data, sizeof(double)*dataLength);
  noiseModel[noiseModelThreadSlot]->NoiseEstimationFFTExpWin(noise_floor, dataLength, NoiseWinSize);
 
  double m, n;
  unsigned long long *ptr;
  for( int i = 0; i < dataLength; i++)
  {
    //Calculate required bit-depth according to the distance to the noise floor
    m = (52 - MIN_MANTISSA_BITS)/( noise_floor[i] * ( NOISE_THRESHOLD_UPPER - NOISE_THRESHOLD_LOWER ) );
    n = MIN_MANTISSA_BITS - m*NOISE_THRESHOLD_LOWER*noise_floor[i];
    resolution_bits = (int)round( m*data[i] + n );
    
    //Limit resolution bits to a double mantissa valid range
    resolution_bits = resolution_bits > 52 ? 52 : resolution_bits;
    resolution_bits = resolution_bits <  1 ?  1 : resolution_bits;
    
    //Apply the mask to double's mantissa
    ptr = (unsigned long long*) (data + i);
    *ptr &= maskLUT_double[resolution_bits-1];  
    
  }

  delete[] noise_floor;
   
}

// [[Rcpp::export]]
void FullImageBitDepthReduction ( StringVector fileNames, 
                                int massChannels, IntegerVector numRows,
                                String dataType, int numOfThreads, 
                                int NoiseWinSize = 16 )
{
  //Copy R data to C arrays
  int numRowsC[fileNames.length()];
  memcpy(numRowsC, numRows.begin(), sizeof(int)*fileNames.length());
  
  MTBitDepthReduction::ImgProcDef myProcParams;
  myProcParams.dataType = dataType;
  myProcParams.fileNames = fileNames;
  myProcParams.massChannels = massChannels;
  myProcParams.numOfThreads = numOfThreads;
  myProcParams.numRows = numRowsC; 
  myProcParams.noiseWinSize = NoiseWinSize;
  
  MTBitDepthReduction myBitReduction(myProcParams);
  myBitReduction.Run();
}

// [[Rcpp::export]]
NumericVector SpectrumBitDepthReduction ( NumericVector data, int NoiseWinSize = 16 )
{
  MTBitDepthReduction::ImgProcDef myProcParams;
  myProcParams.dataType = "double";
  myProcParams.fileNames = StringVector::create();
  myProcParams.massChannels = data.length();
  myProcParams.numOfThreads = 1;
  myProcParams.numRows = new int[1]; 
  myProcParams.noiseWinSize = NoiseWinSize;
  MTBitDepthReduction myBitReduction(myProcParams);
  
  double *x = new double[data.length()];
  memcpy(x, data.begin(), sizeof(double)*data.length());
  myBitReduction.BitDepthReduction(x, data.length(), 1);
  NumericVector res(data.length());
  memcpy(res.begin(), x, sizeof(double)*data.length());
  delete[] x;
  return(res);
}
