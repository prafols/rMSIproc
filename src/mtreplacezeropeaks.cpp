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

#include <Rcpp.h>
#include <cmath>
#include "mtreplacezeropeaks.h"
using namespace Rcpp;

MTReplaceZeroPeaks::MTReplaceZeroPeaks(ImgProcDef imgRunInfo, List peakMatrix):
  ThreadingMsiProc( imgRunInfo.numOfThreads, false, imgRunInfo.basePath, imgRunInfo.fileNames, imgRunInfo.massChannels, imgRunInfo.numRows, imgRunInfo.dataType )
{
  replacedZerosCounters = new int[numOfThreadsDouble];
  for( int i = 0; i < numOfThreadsDouble; i++)
  {
    replacedZerosCounters[i] = 0;
  }
  
  int nrow = as<NumericMatrix>(peakMatrix["intensity"]).rows();
  int ncol = as<NumericMatrix>(peakMatrix["intensity"]).cols();
  pkMatintensity = new NumericMatrix(nrow, ncol);
  pkMatarea = new NumericMatrix(nrow, ncol);
  pkMatsnr = new NumericMatrix(nrow, ncol);
  pkMatmass = new NumericVector(ncol);
  
  *pkMatintensity = as<NumericMatrix>(peakMatrix["intensity"]);
  *pkMatarea = as<NumericMatrix>(peakMatrix["area"]);
  *pkMatsnr = as<NumericMatrix>(peakMatrix["SNR"]);
  *pkMatmass = as<NumericVector>(peakMatrix["mass"]);

  peakObj = new PeakPicking*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    peakObj[i] = new PeakPicking(imgRunInfo.peakWinSize, imgRunInfo.massAxis, imgRunInfo.massChannels, imgRunInfo.peakInterpolationUpSampling );  
  }
  
  //Build the mass index vector
  Rcout<<"Creating the mass index vector...\n";
  mass_index = new int[(*pkMatmass).length()];
  int lastMassId = 1;
  double mass_distance_ant, mass_distance;
  for(int i = 0; i < (*pkMatmass).length(); i++)
  {
    mass_distance_ant = fabs( (*pkMatmass)(i) -  imgRunInfo.massAxis[lastMassId-1]);
    for( int imass = lastMassId; imass <  imgRunInfo.massChannels; imass++)
    {
      mass_distance = fabs( (*pkMatmass)(i) - imgRunInfo.massAxis[imass]);
      if( mass_distance > mass_distance_ant )
      {
        mass_index[i] = imass - 1;
        lastMassId = imass;
        break;
      }
      mass_distance_ant = mass_distance;
    }
  }
}

MTReplaceZeroPeaks::~MTReplaceZeroPeaks()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete peakObj[i];
  }
  delete[] peakObj;
  delete pkMatintensity;
  delete pkMatarea;
  delete pkMatsnr;
  delete pkMatmass;
  delete[] replacedZerosCounters;
  delete[] mass_index;
}

List MTReplaceZeroPeaks::Run()
{
  runMSIProcessingCpp();
  int replacedZerosSum = 0;
  for( int i = 0; i < numOfThreadsDouble; i++)
  {
    replacedZerosSum += replacedZerosCounters[i];
  }
  Rcout << "A total of " << replacedZerosSum << " zeros were replaced\n";
  
  return List::create( Named("mass") = *pkMatmass, Named("intensity") = *pkMatintensity, Named("SNR") = *pkMatsnr, Named("area") = *pkMatarea);
}

void MTReplaceZeroPeaks::ProcessingFunction(int threadSlot)
{
  int cube_row = 0;

  //Look for zeros on current cube corresponding peak matrix rows
  for( int i = CubeNumRows*iCube[threadSlot]; i < (CubeNumRows*iCube[threadSlot] + cubes[threadSlot]->nrows); i++)
  {
    for( int j = 0; j < (*pkMatintensity).ncol(); j ++)
    {
      //Locate the zeros.....
      if( (*pkMatintensity)(i,j) == 0.0 )
      {
        replacedZerosCounters[threadSlot]++;
        //Fill matrix position with proper intensity
        (*pkMatintensity)(i,j) = cubes[threadSlot]->data[cube_row][mass_index[j]];
        (*pkMatarea)(i,j) = peakObj[threadSlot]->predictPeakArea(cubes[threadSlot]->data[cube_row], mass_index[j]);
      }
    }
    cube_row++;
  }
}

// [[Rcpp::export]]
List ReplacePeakMatrixZerosC(List PeakMatrix, String basePath, StringVector fileNames, 
                             NumericVector mass, IntegerVector numRows,
                             String dataType, int numOfThreads, 
                             int WinSize = 10, int InterpolationUpSampling = 10)
{
  //Copy R data to C arrays
  double *massC = new double[mass.length()];
  int numRowsC[fileNames.length()];
  memcpy(massC, mass.begin(), sizeof(double)*mass.length());
  memcpy(numRowsC, numRows.begin(), sizeof(int)*fileNames.length());
  
  MTReplaceZeroPeaks::ImgProcDef imgParams;
  imgParams.basePath = basePath;
  imgParams.dataType = dataType;
  imgParams.fileNames = fileNames;
  imgParams.massAxis = massC;
  imgParams.massChannels = mass.length();
  imgParams.numOfThreads = numOfThreads;
  imgParams.numRows = numRowsC;
  imgParams.peakInterpolationUpSampling = InterpolationUpSampling;
  imgParams.peakWinSize = WinSize;
  
  MTReplaceZeroPeaks pkZeroReplacer(imgParams, PeakMatrix );
  delete[] massC;
  return pkZeroReplacer.Run();
}
