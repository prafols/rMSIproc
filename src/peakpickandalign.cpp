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
#include "peakpickandalign.h"
using namespace Rcpp;

PeakPickAlign::PeakPickAlign():
  ThreadingMsiProc()
{
  //Just init dummy objects to allow the desctructor to work properly
  numOfThreadsDouble = 0;
  numOfPixels = 0;
  peakObj = new PeakPicking*[1];
  alngObj = new LabelFreeAlign*[1];
  mLags  =  new LabelFreeAlign::TLags[1];
  mPeaks =  new PeakPicking::Peaks*[1];
}

PeakPickAlign::PeakPickAlign(ImgProcDef imgRunInfo) : 
  ThreadingMsiProc( imgRunInfo.numOfThreads, imgRunInfo.AlignmentIterations > 0, imgRunInfo.basePath, imgRunInfo.fileNames, imgRunInfo.massChannels, imgRunInfo.numRows, imgRunInfo.dataType )
{
  peakObj = new PeakPicking*[numOfThreadsDouble];
  alngObj = new LabelFreeAlign*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    peakObj[i] = new PeakPicking(imgRunInfo.peakWinSize, imgRunInfo.massAxis, imgRunInfo.massChannels, imgRunInfo.peakInterpolationUpSampling, imgRunInfo.peakSmoothingKernelSize);  
    alngObj[i] = new LabelFreeAlign(imgRunInfo.ref_spectrum, imgRunInfo.massChannels, &fftSharedMutex, imgRunInfo.AlignmentIterations, imgRunInfo.AlignmentMaxShift);
  }
  numberOfCubes =  imgRunInfo.fileNames.length();
  useAlignment = imgRunInfo.AlignmentIterations > 0;
  minSNR = imgRunInfo.SNR;
  numOfPixels = 0;
  binSize = imgRunInfo.tolerance;
  binFilter = imgRunInfo.filter;
  for( int i = 0; i < numberOfCubes; i++)
  {
    numOfPixels +=  imgRunInfo.numRows[i];
  }
  mPeaks =  new PeakPicking::Peaks*[numOfPixels];
  mLags  =  new LabelFreeAlign::TLags[numOfPixels]; 
  bDoBinning = imgRunInfo.performBinning;
}

PeakPickAlign::~PeakPickAlign()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete peakObj[i];
    delete alngObj[i];
  }
  delete[] peakObj;
  delete[] alngObj;

  for( int i = 0; i < numOfPixels; i++ )
  {
    delete mPeaks[i];
  }
  delete[] mPeaks;
  delete[] mLags;
}

List PeakPickAlign::Run()
{
  //Run peak-picking and alignemt in mutli-threading
  runMSIProcessingCpp();

  //Binning
  if(bDoBinning)
  {
    return BinPeaks();
  }
  else
  {
    List pkLst(numOfPixels);
    List pkElem;
    NumericVector pkMass, pkIntensity, pkSNR;
    for(int i = 0; i < numOfPixels; i++)
    {
      Rcout<<"Copying peaks "<< (i+1) <<" of "<< numOfPixels <<"\n";
      pkMass = NumericVector(mPeaks[i]->mass.size());
      pkIntensity = NumericVector(mPeaks[i]->intensity.size());
      pkSNR = NumericVector(mPeaks[i]->SNR.size());
      memcpy(pkMass.begin(), mPeaks[i]->mass.data(), sizeof(double)*mPeaks[i]->mass.size());
      memcpy(pkIntensity.begin(), mPeaks[i]->intensity.data(), sizeof(double)*mPeaks[i]->intensity.size());
      memcpy(pkSNR.begin(), mPeaks[i]->SNR.data(), sizeof(double)*mPeaks[i]->SNR.size());
      pkElem = List::create( Named("mass") = pkMass, Named("intensity") = pkIntensity, Named("SNR") = pkSNR);
      pkLst[i] = pkElem;
    }
    return pkLst;
  }
}


void PeakPickAlign::AppendPeaksAsMatrix(List peaksLst)
{
  //Create a local Peaks to handle all peaks
  int localNumPixels = numOfPixels + (as<NumericMatrix>(peaksLst["intensity"])).rows();
  PeakPicking::Peaks **localPeaks =  new PeakPicking::Peaks*[localNumPixels];
  
  //Copy the current Peaks to the local Peaks object 
  for( int i = 0; i < numOfPixels; i++)
  {
    localPeaks[i] = mPeaks[i];
  }
  delete[] mPeaks; //Delete the old references array, but keeps data pointers in localPeaks
  mPeaks = localPeaks; //Move to localPeaks pointer and set it as current peaks pointer
  
  //Append new peaks 
  for( int i = 0; i < (localNumPixels - numOfPixels); i++)
  {
    mPeaks[ i + numOfPixels] = new PeakPicking::Peaks();
    
    for(int j = 0; j < (as<NumericVector>(peaksLst["mass"])).length(); j++)
    {
      if( (as<NumericMatrix>(peaksLst["intensity"]))(i, j) > 0  )
      {
        mPeaks[ i + numOfPixels]->mass.push_back( (as<NumericVector>(peaksLst["mass"]))(j) ); 
        mPeaks[ i + numOfPixels]->intensity.push_back( (as<NumericMatrix>(peaksLst["intensity"]))(i, j) );
        mPeaks[ i + numOfPixels]->SNR.push_back( (as<NumericMatrix>(peaksLst["SNR"]))(i, j)  );   
        mPeaks[ i + numOfPixels]->area.push_back( (as<NumericMatrix>(peaksLst["area"]))(i, j)  );
      }
    }
  }
  numOfPixels = localNumPixels;
}
  
  
List PeakPickAlign::BinPeaks()
{
  
  //Compute each peak obj TIC, I'll start binning from the highest spectrum which is assumed to have a lower mass error (at least for higher peaks)
  double Tic[numOfPixels];
  for(int i = 0; i < numOfPixels; i++)
  {
    Tic[i] = 0.0;
    for( int j = 0; j < mPeaks[i]->intensity.size(); j++)
    {
      Tic[i] += mPeaks[i]->intensity[j];
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
    
    for( int ip = 0; ip <  mPeaks[iMax]->mass.size(); ip ++)
    {
      binMass.push_back(mPeaks[iMax]->mass[ip]); //Append element to the mass vector (names of bin Matrix)
      binMat.resize(binMat.size() + 1); //Append new column
      binMat[binMat.size() - 1].resize(numOfPixels); //Extenend all new column elements
      binMat[binMat.size() - 1][iMax].intensity = mPeaks[iMax]->intensity[ip]; 
      binMat[binMat.size() - 1][iMax].SNR = mPeaks[iMax]->SNR[ip];
      double numMasses = 1.0; //number of peak masses used to compute each mass centroid
      int countPeaks = 1; //Number of peaks in current column
      
      for( int j = 0; j < numOfPixels; j++)
      {
        if(j != iMax) //Do not process current Peaks obj
        {
          //Find the nearest peak, the nearest peak postion is reatined in iPos var
          double dist, dist_ant = 1e50;
          int iPos = 0;
          for( int imass = 0; imass <  mPeaks[j]->mass.size(); imass++)
          {
            iPos = imass;
            dist = binMass[binMass.size() - 1] - mPeaks[j]->mass[imass];
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
          if(std::abs(dist) <= binSize && mPeaks[j]->mass.size() > 0)
          {
            binMat[binMat.size() - 1][j].intensity = mPeaks[j]->intensity[iPos];
            binMat[binMat.size() - 1][j].SNR = mPeaks[j]->SNR[iPos];
            if(mPeaks[j]->SNR[iPos] >= minSNR)
            {
              //Count peak only if its SNR fits the min SNR value
              countPeaks++;
            }
            
            //Recompute mass centroid using a continuous average
            binMass[binMass.size() - 1] *= numMasses; 
            binMass[binMass.size() - 1] +=  mPeaks[j]->mass[iPos];
            numMasses++;
            binMass[binMass.size() - 1] /= numMasses; 
            
            //Delete datapoint from current peaks
            mPeaks[j]->mass.erase(mPeaks[j]->mass.begin() + iPos);
            mPeaks[j]->intensity.erase(mPeaks[j]->intensity.begin() + iPos);
            mPeaks[j]->SNR.erase(mPeaks[j]->SNR.begin() + iPos);
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
  Rcout<<"Coping data to R object...\n";
  NumericMatrix binMatIntensity(numOfPixels, binMass.size());
  NumericMatrix binMatSNR(numOfPixels, binMass.size());
  NumericMatrix binMatArea(numOfPixels, binMass.size());
  
  //Sort columns by mass
  Rcout<<"Sorting columng by mass...\n";
  NumericVector massSorting(binMass.size());
  memcpy(massSorting.begin(), binMass.begin(), sizeof(double)*binMass.size());
  int sortedInds[binMass.size()];
  double minVal;
  for(int i = 0; i < binMass.size(); i++)
  {
    minVal = std::numeric_limits<double>::max();
    for( int j = 0; j < binMass.size(); j++ )
    {
      if( massSorting[j] < minVal && massSorting[j] > 0 )
      {
        minVal = massSorting[j];
        sortedInds[i] = j;
      }  
    }
    massSorting[ sortedInds[i]  ] = -1; //Mark as sorted
  }
  
  //Copy the mass axis sorting it
  for(int i = 0; i < binMass.size(); i++)
  {
    massSorting[i] = binMass[sortedInds[i]];
  }
  
  //Copy the matrix sorting it
  for( int ir = 0;  ir < numOfPixels; ir++)
  {
    for(int ic = 0; ic < binMass.size(); ic++)
    {
      binMatIntensity(ir, ic ) = binMat[sortedInds[ic]][ir].intensity;
      binMatSNR      (ir, ic ) = binMat[sortedInds[ic]][ir].SNR;
      binMatArea     (ir, ic ) = binMat[sortedInds[ic]][ir].intensity; //TODO currently I'm using area equal intensity for developing but i must calculate area some day, but not here, in peak picking class
    }
  }
  
  //Append lags if alignment was used and zeros otherwise
  NumericVector LagsLow(numOfPixels);
  NumericVector LagsHigh(numOfPixels);
  for( int i = 0; i < numOfPixels; i++)
  {
    if(useAlignment)
    {
      LagsLow[i] = mLags[i].lagLow;
      LagsHigh[i] = mLags[i].lagHigh;
    }
    else
    {
      LagsLow[i] = 0.0;
      LagsHigh[i] = 0.0;
    }
  }
  
  return List::create( Named("mass") = massSorting, Named("intensity") = binMatIntensity, Named("SNR") = binMatSNR, Named("area") = binMatArea, Named("LagLow") = LagsLow, Named("LagHigh") = LagsHigh);
}

void PeakPickAlign::SetBinSize(double value)
{
  binSize = value;
}

void PeakPickAlign::SetBinFilter(double value)
{
  binFilter = value;
}

void PeakPickAlign::ProcessingFunction(int threadSlot)
{
  double snrThreshold = minSNR;
  if(bDoBinning)
  {
    //Using the half of SNR to reatain near-to-min peaks
    snrThreshold *= 0.5;
  }
  
  //Perform alignment and peak-picking of each spectrum in the current loaded cube
  int is = CubeNumRows*iCube[threadSlot];
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    if( useAlignment )
    {
      mLags[is] = alngObj[threadSlot]->AlignSpectrum( cubes[threadSlot]->data[j] );
    }
    mPeaks[is] = peakObj[threadSlot]->peakPicking( cubes[threadSlot]->data[j], snrThreshold ); 
    is++;
  }
}

// [[Rcpp::export]]
List FullImageProcess( String basePath, StringVector fileNames, 
                                NumericVector mass, NumericVector refSpectrum, IntegerVector numRows,
                                String dataType, int numOfThreads, 
                                int AlignmentIterations = 0, int AlignmentMaxShiftPpm = 200,
                                double SNR = 5, int WinSize = 10,
                                int InterpolationUpSampling = 10, int SmoothingKernelSize = 5, 
                                bool doBinning = true, double binningTolerance = 0.05, double binningFilter = 0.9)
{
  
  //Copy R data to C arrays
  double massC[mass.length()];
  double refC[refSpectrum.length()];
  int numRowsC[fileNames.length()];
  memcpy(massC, mass.begin(), sizeof(double)*mass.length());
  memcpy(refC, refSpectrum.begin(), sizeof(double)*refSpectrum.length());
  memcpy(numRowsC, numRows.begin(), sizeof(int)*fileNames.length());
  
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
  myProcParams.AlignmentIterations = AlignmentIterations;
  myProcParams.AlignmentMaxShift = AlignmentMaxShiftPpm;
  myProcParams.SNR = SNR;
  myProcParams.tolerance = binningTolerance;
  myProcParams.filter = binningFilter;
  myProcParams.ref_spectrum = refC;
  myProcParams.performBinning = doBinning;
  
  PeakPickAlign myPeakPicking(myProcParams);
  return myPeakPicking.Run();
}

// [[Rcpp::export]]
List MergePeakMatricesC( List PeakMatrices, double binningTolerance = 0.05, double binningFilter = 0.01 )
{
  PeakPickAlign myPeakPicking;
  for( int i = 0; i < PeakMatrices.length(); i++)
  {
    Rcout<<"Merging peak matrix "<<(i+1)<<" of "<<PeakMatrices.length()<<"\n";
    myPeakPicking.AppendPeaksAsMatrix(PeakMatrices[i]);  
  }
  myPeakPicking.SetBinFilter(binningFilter);
  myPeakPicking.SetBinSize(binningTolerance);
  return myPeakPicking.BinPeaks();
}
