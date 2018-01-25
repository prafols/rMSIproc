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

#ifndef LABEL_FREE_ALIGN_H
  #define LABEL_FREE_ALIGN_H

#include <Rcpp.h>
#include <fftw3.h>
#include <boost/thread.hpp>

#define BOTTOM_SPECTRUM 0
#define CENTER_SPECTRUM 1
#define TOP_SPECTRUM 2

class LabelFreeAlign
{
  public:
    //spectraSplit the low/high part of spectra to keep (the resting points to unit will be removed).
    LabelFreeAlign(double *mass, double *ref_spectrum, int numOfPoints, bool bilinear, boost::mutex *sharedMutex, int iterations = 3, 
                   double lagRefLow = 0.1, double lagRefMid = 0.5, double lagRefHigh = 0.9,
                   double lagLimitppm = 200, int fftOverSampling = 2, double winSizeRelative = 0.6);
    ~LabelFreeAlign();
    
    //Data accessors to internal vars to test this class
    Rcpp::NumericVector getHannWindow();
    Rcpp::NumericVector getHannWindowCenter();
    Rcpp::NumericVector getRefLowFFT();
    Rcpp::NumericVector getRefCenterFFT();
    Rcpp::NumericVector getRefHighFFT();
    
    typedef struct
    {
      double lagHigh;
      double lagMid;
      double lagLow;
    }TLags;
    
    //Algin a given spectrum to the reference specified in the constructor. 
    //The used High and Low lags are returned as a TLags structure.
    TLags AlignSpectrum(double *data );

  private:
    void ComputeRef(double *data_ref, int spectrumPart);
    void ZeroPadding(double *data,bool reverse, int targetSize, int dataSize);
    void CopyData2Window(double *data_int, double *data_out,  int spectrumPart);
    void TimeWindow(double *data, int spectrumPart);
    int FourierBestCor(double *data, double *ref);
    void FourierLinerScaleShift(double *data, double scaling, double shift);
  
    int dataLength; //Number of points used in each spectrum
    int WinLength; //Number of points of spectrum retained in hanning window
    int FFT_Size; //Number of points used for fft
    int FFT_Size_SH; //Number of points used for the first fft of time scale/shift
    fftw_plan fft_pdirect;
    fftw_plan fft_pinvers;
    fftw_plan fft_pdshiftScale;
    double *fft_direct_in;
    double *fft_direct_out;
    double *fft_inverse_in;
    double *fft_inverse_out;
    fftw_complex *fft_direct_shiftScale_in;
    fftw_complex *fft_direct_shiftScale_out;

    //Mem space to store pre-computed reference FFT space values
    double *fft_ref_low;
    double *fft_ref_center;
    double *fft_ref_high;
    double *HannWindow;
    double *HannWindowCenter;
    double lagMaxLow;
    double lagMaxMid;
    double lagMaxHigh;
    bool bBilinear;
    int AlignIterations;
    double RefLagLow;
    double RefLagMid;
    double RefLagHigh;
    int FFTScaleShiftOverSampling;
    
    boost::mutex *fftwMtx; //Lock mechanism for signalling the FourierLinerScaleShift fftw non-thread-safe calls
};

#endif