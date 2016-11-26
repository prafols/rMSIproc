/*************************************************************************
 * This file is part of rMSIproc.
 *
 * rMSIproc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * rMSIproc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rMSIproc.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "peakpickandalign.h"
using namespace Rcpp;

PeakPickAlign::PeakPickAlign(ImgProcDef imgRunInfo) : 
  ThreadingMsiProc( imgRunInfo.numOfThreads, imgRunInfo.runAlignment, imgRunInfo.basePath, imgRunInfo.fileNames, imgRunInfo.massChannels, imgRunInfo.numRows, imgRunInfo.dataType )
{
  peakObj = new PeakPicking(imgRunInfo.peakWinSize, imgRunInfo.massAxis, imgRunInfo.massChannels, imgRunInfo.peakInterpolationUpSampling, imgRunInfo.peakSmoothingKernelSize);
  numberOfCubes =  imgRunInfo.fileNames.length();
  useAlignment = imgRunInfo.runAlignment;
  minSNR = imgRunInfo.SNR;
  numOfPixels = 0;
  binSize = imgRunInfo.tolerance;
  binFilter = imgRunInfo.filter;
  for( int i = 0; i < numberOfCubes; i++)
  {
    numOfPixels +=  imgRunInfo.numRows[i];
  }
}

PeakPickAlign::~PeakPickAlign()
{
  delete peakObj;
}

List PeakPickAlign::Run()
{
  CrMSIDataIO::DataCube *dc;
  PeakPicking::Peaks *mPeaks[numOfPixels];
  int is = 0; //The spectrum indexer
  
  //TODO this is a mono-thread implementation
  for( int i = 0; i < numberOfCubes; i++)
  {
    Rcout<<"Processing cube "<<(i + 1)<<" of "<<numberOfCubes<<"\n";
    
    //Load a datacube
    dc = ioObj->loadDataCube(i);
    
    if( useAlignment )
    {
      //TODO call alignment object here
      //Alignment will overwirte the ramdisk so WARNING! use te algiend spectra to perform peakpicking!!!!!
    }
    
    //Perform peak-picking of each spectrum in the current loaded cube
    for( int j = 0; j < dc->nrows; j++)
    {
      mPeaks[is] = peakObj->peakPicking( dc->data[j], minSNR );
      is++;
    }
    
    //Clear memory of current cube
    ioObj->freeDataCube(dc);
  }
  
  //Binning 
  List mBins = BinPeaks(mPeaks);
  
  //Free memory
  for( int i = 0; i < numOfPixels; i++ )
  {
    delete mPeaks[i];
  }
  Rcout<<"Processing complete\n";
  
  return mBins;
}

List PeakPickAlign::BinPeaks(PeakPicking::Peaks **peaks)
{
  
  //Compute each peak obj TIC, I'll start binning from the highest spectrum which is assumed to have a lower mass error (at least for higher peaks)
  double Tic[numOfPixels];
  for(int i = 0; i < numOfPixels; i++)
  {
    Tic[i] = 0.0;
    for( int j = 0; j < peaks[i]->intensity.size(); j++)
    {
      Tic[i] += peaks[i]->intensity[j];
    }
  }
  
  //Apply binning starting with most intens spectra
  double dMax;
  int iMax;
  std::vector<std::vector<TBin> > binMat; //A matrix containing binning values
  NumericVector binMass; //The mass name for each matrix column
  Rcout<<"Binning...\n";
  int iCount = 0;
  while(true)
  {
    Rcout<<iCount<<"/"<<numOfPixels<<"\n";
    iCount++;
    //Locate current most intens spectrum
    iMax = -1;
    dMax = -1.0;
    for(int i = 0; i < numOfPixels; i++)
    {
      if(Tic[i] > dMax && Tic[i] >= 0)
      {
        dMax = Tic[i];
        iMax = i;
      }
    }
    
    //End condition
    if(iMax == -1)
    {
      break;
    }
    
    Tic[iMax] = -1.0; //Mark current spectrum as processed by deleting its TIC
    
    for( int ip = 0; ip <  peaks[iMax]->mass.size(); ip ++)
    {
      binMass.push_back(peaks[iMax]->mass[ip]); //Append element to the mass vector (names of bin Matrix)
      binMat.resize(binMat.size() + 1); //Append new column
      binMat[binMat.size() - 1].resize(numOfPixels); //Extenend all new column elements
      binMat[binMat.size() - 1][iMax].intensity = peaks[iMax]->intensity[ip]; 
      binMat[binMat.size() - 1][iMax].SNR = peaks[iMax]->SNR[ip];
      double numMasses = 1.0; //number of peak masses used to compute each mass centroid
      int countPeaks = 1; //Number of peaks in current column
      
      for( int j = 0; j < numOfPixels; j++)
      {
        if(j != iMax) //Do not process current Peaks obj
        {
          //Find the nearest peak, the nearest peak postion is reatined in iPos var
          double dist, dist_ant = 1e50;
          int iPos = 0;
          for( int imass = 0; imass <  peaks[j]->mass.size(); imass++)
          {
            iPos = imass;
            dist = binMass[binMass.size() - 1] - peaks[j]->mass[imass];
            if(dist <= 0)
            {
              if(std::abs(dist_ant) < std::abs(dist))
              {
                iPos = imass - 1;
                dist = dist_ant;
              }
              break; //avoid all the loop
            }
            dist_ant = dist;
          }
          
          //Check if is in the same mass bin and fill the matrix value accordingly
          if(std::abs(dist) <= binSize && peaks[j]->mass.size() > 0)
          {
            binMat[binMat.size() - 1][j].intensity = peaks[j]->intensity[iPos];
            binMat[binMat.size() - 1][j].SNR = peaks[j]->SNR[iPos];
            countPeaks++;
            
            //Recompute mass centroid using a continuous average
            binMass[binMass.size() - 1] *= numMasses; 
            binMass[binMass.size() - 1] +=  peaks[j]->mass[iPos];
            numMasses++;
            binMass[binMass.size() - 1] /= numMasses; 
            
            //Delete datapoint from current peaks
            peaks[j]->mass.erase(peaks[j]->mass.begin() + iPos);
            peaks[j]->intensity.erase(peaks[j]->intensity.begin() + iPos);
            peaks[j]->SNR.erase(peaks[j]->SNR.begin() + iPos);
          }
          else
          {
            binMat[binMat.size() - 1][j].intensity = 0.0;
            binMat[binMat.size() - 1][j].SNR = 0.0;
          }
        }
      }
      
      //Delete last column if it does not fit minimum filter criterion
      if(countPeaks < (int)(numOfPixels * binFilter))
      {
        binMass.erase(binMass.begin() + binMass.size() - 1);
        binMat.erase(binMat.begin() + binMat.size() - 1);
      }
      
    }
  }
  Rcout<<"Bining complete with a total number of "<<binMass.size()<<" bins\n";
  
  //Copy data to R matrices
  Rcout<<"Coping data to R object...";
  NumericMatrix binMatIntensity(numOfPixels, binMass.size());
  NumericMatrix binMatSNR(numOfPixels, binMass.size());
  for( int ir = 0;  ir < numOfPixels; ir++)
  {
    for(int ic = 0; ic < binMass.size(); ic++)
    {
      binMatIntensity(ir, ic) = binMat[ic][ir].intensity;
      binMatSNR(ir, ic) = binMat[ic][ir].SNR;
    }
  }
  return List::create( Named("mass") = binMass, Named("intensity") = binMatIntensity, Named("SNR") = binMatSNR);
}

// [[Rcpp::export]]
List FullImageProcess( String basePath, StringVector fileNames, 
                                NumericVector mass, IntegerVector numRows, 
                                String dataType, int numOfThreads, 
                                bool runAlignment = false, double SNR = 5, int WinSize = 10,
                                int InterpolationUpSampling = 10, int SmoothingKernelSize = 5, 
                                double binningTolerance = 0.05, double binningFilter = 0.9)
{
  
  //Copy R data to C arrays
  double massC[mass.length()];
  int numRowsC[fileNames.length()];
  std::memcpy(massC, mass.begin(), sizeof(double)*mass.length());
  std::memcpy(numRowsC, numRows.begin(), sizeof(int)*fileNames.length());
  
  PeakPickAlign::ImgProcDef myProcParams;
  myProcParams.basePath = basePath;
  myProcParams.dataType = dataType;
  myProcParams.fileNames = fileNames;
  myProcParams.massAxis = massC;
  myProcParams.massChannels = mass.length();
  myProcParams.numOfThreads = numOfThreads;
  myProcParams.numRows = numRowsC; 
  myProcParams.peakInterpolationUpSampling  = InterpolationUpSampling;
  myProcParams.peakSmoothingKernelSize = SmoothingKernelSize;
  myProcParams.peakWinSize = WinSize;
  myProcParams.runAlignment = runAlignment;
  myProcParams.SNR = SNR;
  myProcParams.tolerance = binningTolerance;
  myProcParams.filter = binningFilter;
  
  PeakPickAlign myPeakPicking(myProcParams);
  return myPeakPicking.Run();
}
