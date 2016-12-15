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

#ifndef PEAKPICK_AND_ALIGN_H
  #define PEAKPICK_AND_ALIGN_H
#include <Rcpp.h>
#include <boost/thread.hpp>
#include "peakpicking.h"
#include "labelfreealign.h"
#include "threadingmsiproc.h"


class PeakPickAlign : public ThreadingMsiProc 
{
  public:
    //Data structur used to completly define the processing pipeline
    typedef struct
    {
      int peakWinSize; //Windows size used for peak-picking
      double *massAxis; //Array containing the mass axis
      double *ref_spectrum; //Reference spectrum used for alignment
      int massChannels; //Number of points in the mass axis array
      int peakInterpolationUpSampling; //Upsampling value for peak interpolation
      int peakSmoothingKernelSize; //SavitzkyGolay kernel size for smoothing
      double SNR; //Minimum peak signal to noise ratio
      Rcpp::String basePath; //Full path where ramdisks are stored
      Rcpp::StringVector fileNames; //Filname of each ramdisk file
      int *numRows; //An array containing the number of rows stored in each ramdisk file. The length ot this array is the length of fileNames
      Rcpp::String dataType; //A string with the data type
      int AlignmentIterations;
      int numOfThreads;
      double tolerance;
      double filter;
      bool performBinning;
    }ImgProcDef;  
    
    PeakPickAlign();
    PeakPickAlign(ImgProcDef imgRunInfo);
    ~PeakPickAlign();
    
    //Exectur a full imatge processing using threaded methods
    Rcpp::List Run(); 
    
    //Appends a list of peaks to the current list of peaks (mPeaks)
    //Peaks are input as a binned matrix
    void AppendPeaksAsMatrix(Rcpp::List peaksLst);
    
    //Perfomr peak binning over mPeaks, this is mono-thread implemented.
    Rcpp::List BinPeaks();
    void SetBinSize(double value);
    void SetBinFilter(double value);
    
  private:
    PeakPicking **peakObj;
    LabelFreeAlign **alngObj;
    int numberOfCubes;
    bool useAlignment;
    int numOfPixels;
    double minSNR;
    double binSize;
    double binFilter;
    PeakPicking::Peaks **mPeaks; //A place to store peaks objects outside threaded space
    LabelFreeAlign::TLags *mLags; //A place to store alignment lags
    bool bDoBinning;
    
    typedef struct
    {
      double intensity;
      double SNR;
    }TBin;
    
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
    
    boost::mutex fftSharedMutex;
};
  
  
#endif