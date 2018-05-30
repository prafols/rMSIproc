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
#include <vector>
#include <limits>
#include "mtpeakpicking.h"
#include "peakbinning.h"
using namespace Rcpp;

MTPeakPicking::MTPeakPicking(ImgProcDef imgRunInfo) : 
  ThreadingMsiProc( imgRunInfo.numOfThreads, false, imgRunInfo.fileNames, imgRunInfo.massChannels, imgRunInfo.numRows, imgRunInfo.dataType )
{
  peakObj = new PeakPicking*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    peakObj[i] = new PeakPicking(imgRunInfo.peakWinSize, imgRunInfo.massAxis, imgRunInfo.massChannels, imgRunInfo.peakInterpolationUpSampling );  
  }
  
  minSNR = imgRunInfo.SNR;
  numOfPixels = 0;
  for( int i = 0; i < imgRunInfo.fileNames.length(); i++)
  {
    numOfPixels +=  imgRunInfo.numRows[i];
  }
  mPeaks =  new PeakPicking::Peaks*[numOfPixels];
}

MTPeakPicking::~MTPeakPicking()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete peakObj[i];
  }
  delete[] peakObj;

  for( int i = 0; i < numOfPixels; i++ )
  {
    delete mPeaks[i];
  }
  delete[] mPeaks;
}

void MTPeakPicking::Run()
{
  //Run peak-picking in mutli-threading
  runMSIProcessingCpp();
}

Rcpp::List MTPeakPicking::ExportPeakList()
{
  List pkLst(numOfPixels);
  NumericVector pkMass, pkIntensity, pkSNR, pkArea, pkBinSize;
  for(int i = 0; i < numOfPixels; i++)
  {
    Rcout<<"Copying peaks "<< (i+1) <<" of "<< numOfPixels <<"\n";
    pkMass = NumericVector(mPeaks[i]->mass.size());
    pkIntensity = NumericVector(mPeaks[i]->intensity.size());
    pkSNR = NumericVector(mPeaks[i]->SNR.size());
    pkArea = NumericVector(mPeaks[i]->area.size());
    pkBinSize = NumericVector(mPeaks[i]->binSize.size());
    memcpy(pkMass.begin(), mPeaks[i]->mass.data(), sizeof(double)*mPeaks[i]->mass.size());
    memcpy(pkIntensity.begin(), mPeaks[i]->intensity.data(), sizeof(double)*mPeaks[i]->intensity.size());
    memcpy(pkSNR.begin(), mPeaks[i]->SNR.data(), sizeof(double)*mPeaks[i]->SNR.size());
    memcpy(pkArea.begin(), mPeaks[i]->area.data(), sizeof(double)*mPeaks[i]->area.size());
    memcpy(pkBinSize.begin(), mPeaks[i]->binSize.data(), sizeof(double)*mPeaks[i]->binSize.size());
    pkLst[i] = List::create( Named("mass") = pkMass, Named("intensity") = pkIntensity, Named("SNR") = pkSNR, Named("area") = pkArea, Named("binSize") = pkBinSize);
  }
  return pkLst;
}

void MTPeakPicking::ProcessingFunction(int threadSlot)
{
  //Perform peak-picking of each spectrum in the current loaded cube
  int is = CubeFirstRowID[iCube[threadSlot]];
  //The same will happen for Alignment
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    mPeaks[is] = peakObj[threadSlot]->peakPicking( cubes[threadSlot]->data[j], minSNR ); 
    is++;
  }
}

PeakPicking::Peaks **MTPeakPicking::getPeakList()
{
  return mPeaks;
}

// [[Rcpp::export]]
List FullImagePeakPicking ( StringVector fileNames, 
                                NumericVector mass, IntegerVector numRows,
                                String dataType, int numOfThreads, 
                                double SNR = 5, int WinSize = 10, int InterpolationUpSampling = 10, 
                                bool doBinning = true, double binningTolerance = 100, double binningFilter = 0.9, bool binningIn_ppm = true,
                                bool exportPeakList = false)
{
  //Copy R data to C arrays
  double *massC = new double[mass.length()];
  int numRowsC[fileNames.length()];
  memcpy(massC, mass.begin(), sizeof(double)*mass.length());
  memcpy(numRowsC, numRows.begin(), sizeof(int)*fileNames.length());
  
  MTPeakPicking::ImgProcDef myProcParams;
  myProcParams.dataType = dataType;
  myProcParams.fileNames = fileNames;
  myProcParams.massAxis = massC;
  myProcParams.massChannels = mass.length();
  myProcParams.numOfThreads = numOfThreads;
  myProcParams.numRows = numRowsC; 
  myProcParams.peakInterpolationUpSampling  = InterpolationUpSampling;
  myProcParams.peakWinSize = WinSize;
  myProcParams.SNR = SNR;
  
  MTPeakPicking myPeakPicking(myProcParams);
  delete[] massC;
  
  myPeakPicking.Run();
  PeakPicking::Peaks **CpeakList_ptr = myPeakPicking.getPeakList();
    
  //Call binning if it is enabled
  List peakMat;
  if(doBinning)
  {
    int pixelCount = 0;
    for( int i = 0; i < fileNames.length(); i++)
    {
      pixelCount +=  numRowsC[i];
    }
    PeakBinning myPeakBinning(CpeakList_ptr, pixelCount, binningTolerance, binningIn_ppm, binningFilter);
    peakMat = myPeakBinning.BinPeaks(); 
  }

  //Export the peak list without binning for imzML centroid convertion
  List peakList;
  if( exportPeakList )
  {
    peakList = myPeakPicking.ExportPeakList();
  }
  
  //The peak list object is automatically destroyed by the MTPeakPicking destructor
  return List::create( Named("PeakMatrix") = peakMat, Named("PeakList") = peakList);
}
