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

class LabelFreeAlign
{
  public:
    //spectraSplit the low/high part of spectra to keep (the resting points to unit will be removed).
    LabelFreeAlign(double *ref_spectrum, int numOfPoints,  boost::mutex *sharedMutex, int iterations = 3, double spectraSplit = 0.6, double lagLimitppm = 200);
    ~LabelFreeAlign();
    
    //Data accessors to internal vars to test this class
    Rcpp::NumericVector getHannWindow();
    Rcpp::NumericVector getRefLowFFT();
    Rcpp::NumericVector getRefHighFFT();
    
    typedef struct
    {
      double lagHigh;
      double lagLow;
    }TLags;
    
    //Algin a given spectrum to the reference specified in the constructor. 
    //The used High and Low lags are returned as a TLags structure.
    TLags AlignSpectrum(double *data );

 private:
    void ComputeRef(double *data_ref, bool bHigh);
    void ZeroPadding(double *data,bool reverse, int targetSize, int dataSize);
    void FourierLinerScaleShift(double *data, double scaling, double shift);
    
    //bLow if is true then apply the hann win to the low part of spectrum
    void TimeWindow(double *data, bool bHigh);
    int FourierBestCor(double *data, double *ref);
  
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
    double *fft_ref_high;
    double *HannWindow;
    double lagMax;
    int AlignIterations;
    
    boost::mutex *fftwMtx; //Lock mechanism for signalling the FourierLinerScaleShift fftw non-thread-safe calls
};

#endif