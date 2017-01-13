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
#include "mtlabelfreealign.h"
using namespace Rcpp;

MTLabelFreeAlign::MTLabelFreeAlign(ImgProcDef imgRunInfo) : 
  ThreadingMsiProc( imgRunInfo.numOfThreads, true, imgRunInfo.basePath, imgRunInfo.fileNames, imgRunInfo.massChannels, imgRunInfo.numRows, imgRunInfo.dataType )
{
  alngObj = new LabelFreeAlign*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    alngObj[i] = new LabelFreeAlign(imgRunInfo.ref_spectrum, imgRunInfo.massChannels, &fftSharedMutex, imgRunInfo.AlignmentIterations, imgRunInfo.AlignmentMaxShift);
  }
  numOfPixels = 0;
  for( int i = 0; i < imgRunInfo.fileNames.length(); i++)
  {
    numOfPixels +=  imgRunInfo.numRows[i];
  }
  mLags  =  new LabelFreeAlign::TLags[numOfPixels]; 
}

MTLabelFreeAlign::~MTLabelFreeAlign()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete alngObj[i];
  }
  delete[] alngObj;
  delete[] mLags;
}

List MTLabelFreeAlign::Run()
{
  //Run peak-picking and alignemt in mutli-threading
  runMSIProcessingCpp();

  //Retrun first iteration lags
  NumericVector LagsLow(numOfPixels);
  NumericVector LagsHigh(numOfPixels);
  for( int i = 0; i < numOfPixels; i++)
  {
    LagsLow[i] = mLags[i].lagLow;
    LagsHigh[i] = mLags[i].lagHigh;
  }
  
  return List::create( Named("LagLow") = LagsLow, Named("LagHigh") = LagsHigh);
}

void MTLabelFreeAlign::ProcessingFunction(int threadSlot)
{
  //Perform alignment of each spectrum in the current loaded cube
  int is = CubeNumRows*iCube[threadSlot];
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    mLags[is] = alngObj[threadSlot]->AlignSpectrum( cubes[threadSlot]->data[j] );
    is++;
  }
}

// [[Rcpp::export]]
List FullImageAlign( String basePath, StringVector fileNames, 
                                NumericVector refSpectrum, IntegerVector numRows,
                                String dataType, int numOfThreads, 
                                int AlignmentIterations = 3, int AlignmentMaxShiftPpm = 200 )
{
  
  //Copy R data to C arrays
  double refC[refSpectrum.length()];
  int numRowsC[fileNames.length()];
  memcpy(refC, refSpectrum.begin(), sizeof(double)*refSpectrum.length());
  memcpy(numRowsC, numRows.begin(), sizeof(int)*fileNames.length());
  
  MTLabelFreeAlign::ImgProcDef myProcParams;
  myProcParams.basePath = basePath;
  myProcParams.dataType = dataType;
  myProcParams.fileNames = fileNames;
  myProcParams.massChannels = refSpectrum.length();
  myProcParams.numOfThreads = numOfThreads;
  myProcParams.numRows = numRowsC; 
  myProcParams.AlignmentIterations = AlignmentIterations;
  myProcParams.AlignmentMaxShift = AlignmentMaxShiftPpm;
  myProcParams.ref_spectrum = refC;
  
  MTLabelFreeAlign myAlignment(myProcParams);
  return myAlignment.Run();
}
