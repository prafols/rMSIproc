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
 
#include "labelfreealign.h"

using namespace Rcpp;

LabelFreeAlign::LabelFreeAlign(double *ref_spectrum, int numOfPoints,  boost::mutex *sharedMutex, double spectraSplit,double lagLimitppm)
{
  dataLength = numOfPoints;
  WinLength = (int)round(spectraSplit * (double)dataLength);
  FFT_Size_SH = (int)pow(2.0, std::ceil(log2(dataLength)));
  
  //Pre-compute the Hanning Windows
  HannWindow = new double[WinLength];
  for(int i = 0; i < WinLength; i++)
  {
    HannWindow[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i + 1.0) / ( 2.0 * WinLength ) ));
  }
  
  FFT_Size = (int)pow(2.0, std::ceil(log2(WinLength)));
  fft_direct_in = fftw_alloc_real(FFT_Size);
  fft_direct_out = fftw_alloc_real(FFT_Size);
  fft_inverse_in = fftw_alloc_real(FFT_Size);
  fft_inverse_out = fftw_alloc_real(FFT_Size);
  fft_direct_shiftScale_in = fftw_alloc_complex(FFT_Size_SH);
  fft_direct_shiftScale_out = fftw_alloc_complex(FFT_Size_SH);
  
  fft_pdirect = fftw_plan_r2r_1d(FFT_Size, fft_direct_in, fft_direct_out, FFTW_R2HC, FFTW_ESTIMATE);
  fft_pinvers = fftw_plan_r2r_1d(FFT_Size, fft_inverse_in, fft_inverse_out, FFTW_HC2R, FFTW_ESTIMATE);
  fft_pdshiftScale = fftw_plan_dft_1d(FFT_Size_SH, fft_direct_shiftScale_in, fft_direct_shiftScale_out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  fft_ref_low = new  double[FFT_Size];
  fft_ref_high = new  double[FFT_Size];
  ComputeRef(ref_spectrum, true);
  ComputeRef(ref_spectrum, false);
  
  //Lag limits relative to number of points and FFT siexe
  lagMax = (lagLimitppm/1.0e6)*(double)FFT_Size;
  
  //The mutext shared for all threads
  fftwMtx = sharedMutex;
}

LabelFreeAlign::~LabelFreeAlign()
{
  fftw_destroy_plan(fft_pdirect);
  fftw_destroy_plan(fft_pinvers);
  fftw_destroy_plan(fft_pdshiftScale);
  fftw_free(fft_direct_in); 
  fftw_free(fft_direct_out);
  fftw_free(fft_inverse_in); 
  fftw_free(fft_inverse_out);
  fftw_free(fft_direct_shiftScale_in);
  fftw_free(fft_direct_shiftScale_out);
  delete[] fft_ref_low;
  delete[] fft_ref_high;
  delete[] HannWindow;
}

void LabelFreeAlign::ComputeRef(double *data_ref, bool bHigh)
{
  if(bHigh)
  {
    memcpy( (fft_direct_in + (FFT_Size - WinLength)), data_ref + (dataLength - WinLength), sizeof(double)*WinLength);
  }
  else
  {
    memcpy(fft_direct_in, data_ref, sizeof(double)*WinLength);
  }
  TimeWindow(fft_direct_in,  bHigh);
  ZeroPadding(fft_direct_in, bHigh, FFT_Size, WinLength);
  fftw_execute(fft_pdirect);
  double *data_ptr;
  if(bHigh)
  {
    data_ptr = fft_ref_high;
  }
  else
  {
    data_ptr = fft_ref_low;
  }
  memcpy(data_ptr, fft_direct_out, sizeof(double)*FFT_Size);
  
  //Compute Conj
  for( int i = ((FFT_Size + 1)/2 - 1) ; i < FFT_Size; i++)
  {
    data_ptr[i] *= (-1.0); 
  }
}

void LabelFreeAlign::ZeroPadding(double *data,bool reverse, int targetSize, int dataSize)
{
  if(dataSize == targetSize)
  {
    //No padding needed
    return;
  }
  if(reverse)
  {
    for( int i = 0; i < (targetSize - dataSize); i++)
    {
      data[i] = 0.0;
    }
  }
  else
  {
    for( int i = dataSize; i < targetSize; i++)
    {
      data[i] = 0.0;
    }
  }
}

void LabelFreeAlign::TimeWindow(double *data, bool bHigh)
{
  for( int i = 0; i < WinLength; i++)
  {
    if(bHigh)
    {
      data[i + (FFT_Size - WinLength)] *= HannWindow[i]; 
    }
    else
    {
      data[i] *= HannWindow[WinLength - i - 1]; 
    }
  }
}

LabelFreeAlign::TLags LabelFreeAlign::AlignSpectrum(double *data )
{
  //Hanning Windowing
  double topWin_data[FFT_Size];
  double botWin_data[FFT_Size];
  memcpy( (topWin_data + (FFT_Size - WinLength)), data + (dataLength - WinLength), sizeof(double)*WinLength);
  memcpy(botWin_data, data, sizeof(double)*WinLength);
  TimeWindow(topWin_data, true);
  TimeWindow(botWin_data, false);
  
  //Zero-padding 2 improve fft performance
  ZeroPadding(topWin_data, true, FFT_Size, WinLength);
  ZeroPadding(botWin_data, false, FFT_Size, WinLength);
  
  //Get lags
  TLags lags;
  lags.lagLow = (double)FourierBestCor(botWin_data, fft_ref_low);
  lags.lagHigh = (double)FourierBestCor(topWin_data, fft_ref_high);
  
  //Limit lags
  lags.lagLow = abs(lags.lagLow) > lagMax ? 0.0 : lags.lagLow;
  lags.lagHigh = abs(lags.lagHigh) > lagMax ? 0.0 : lags.lagHigh; 

  //Aligment constansts mz'(n) = K*mz(n) + Sh
  const double Rl = 0.0; //Currently I only use zero as the first reference
  const double Rh =  0.9 * (double)FFT_Size; //Currently I only use N as the first reference
  const double K = (Rh + lags.lagHigh - Rl - lags.lagLow)/(Rh - Rl); //New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
  const double Sh = (Rh* lags.lagLow - Rl*lags.lagHigh)/(Rh - Rl); //If scaling is performed before shift. New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
  
  //Apply the scaling and shift to original data, after that the pointer to data contains the aligned sptectrum
  FourierLinerScaleShift(data, K, Sh);
    
  return lags;
}

void LabelFreeAlign::FourierLinerScaleShift(double *data, double scaling, double shift)
{
  const int NewFFT_Size = FFT_Size_SH + (int)round( (scaling - 1.0) * (double)FFT_Size_SH);
  fftwMtx->lock();
  fftw_complex *fft_odd_in = fftw_alloc_complex(NewFFT_Size);
  fftw_complex *fft_odd_out = fftw_alloc_complex(NewFFT_Size);
  fftw_plan fft_pOddInvers = fftw_plan_dft_1d(NewFFT_Size, fft_odd_in, fft_odd_out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwMtx->unlock();
  
  //Copy data as complex number and add padding zeros and calculate spectrum Max
  double maxInit = 0.0;
  for(int i = 0; i < FFT_Size_SH; i++)
  {
    if(i < dataLength)
    {
      fft_direct_shiftScale_in[i][0] = data[i];
    }
    else
    {
      fft_direct_shiftScale_in[i][0] = 0.0;
    }
    fft_direct_shiftScale_in[i][1] = 0.0; //The imaginary part is always zero because input es real
    
    maxInit = fft_direct_shiftScale_in[i][0] > maxInit? fft_direct_shiftScale_in[i][0] : maxInit;
  }
  fftw_execute(fft_pdshiftScale);
  
  //Copy data to fft_odd_in and apply scaling and shift
  double arg, Re, Im, auxRe, auxIm;
  for( int i = 0; i < NewFFT_Size; i++ )
  {
    //Scaling in FFT domain
    if( i < FFT_Size_SH/2 && i < NewFFT_Size/2)
    {
      fft_odd_in[i][0] = fft_direct_shiftScale_out[i][0]; //Copy real part
      fft_odd_in[i][1] = fft_direct_shiftScale_out[i][1]; //Copy imaginary part
    }
    else if( i >= (NewFFT_Size - FFT_Size_SH/2))
    {
      fft_odd_in[i][0] = fft_direct_shiftScale_out[i - NewFFT_Size + FFT_Size_SH][0]; //Copy real part
      fft_odd_in[i][1] = fft_direct_shiftScale_out[i - NewFFT_Size + FFT_Size_SH][1]; //Copy imaginary part
    }
    else
    {
      //Zero padding
      fft_odd_in[i][0] = 0.0;
      fft_odd_in[i][1] = 0.0;
    }
    
    //Shift in FFT domain
    arg = -2.0*M_PI*(1.0/((double)NewFFT_Size))*((double)i + 1.0)*shift;
    Re = cos(arg);
    Im = sin(arg);
    auxRe = fft_odd_in[i][0] * Re - fft_odd_in[i][1]*Im;
    auxIm = fft_odd_in[i][0] * Im + fft_odd_in[i][1]*Re;
    fft_odd_in[i][0] = auxRe;
    fft_odd_in[i][1] = auxIm;
  }
  fftw_execute(fft_pOddInvers);

  //Copy data to output vector 
  double maxEnd = 0.0;
  for(int i = 0; i < dataLength; i++)
  {
    if( i < NewFFT_Size )
    {
      data[i] = sqrt(fft_odd_out[i][0]*fft_odd_out[i][0] + fft_odd_out[i][1]*fft_odd_out[i][1]); //Module
    }
    else
    {
      data[i] = 0.0;
    }
    maxEnd = data[i] > maxEnd? data[i] : maxEnd;
  }
  
  //Compensate gain
  for(int i = 0; i < dataLength; i++)
  {
    data[i] /= maxEnd/maxInit;
  }
  
  fftwMtx->lock();
  fftw_destroy_plan(fft_pOddInvers);
  fftw_free(fft_odd_in);
  fftw_free(fft_odd_out);
  fftwMtx->unlock();
}

int LabelFreeAlign::FourierBestCor(double *data, double *ref)
{
  memcpy(fft_direct_in, data, sizeof(double)*FFT_Size);
  fftw_execute(fft_pdirect);
  
  //Mult fft complex values, the ref is assumed already Conj (this is automatically done by ComputeRef method)
  for( int i = 0; i <= FFT_Size/2; i++)
  {
    if( i > 0 && i < FFT_Size/2)
    {
      fft_inverse_in[i] = ref[i] * fft_direct_out[i] - ref[FFT_Size - i] * fft_direct_out[FFT_Size - i];
      fft_inverse_in[FFT_Size - i] = ref[i] * fft_direct_out[FFT_Size - i] + ref[FFT_Size - i] * fft_direct_out[i];
    }
    else
    {
      fft_inverse_in[i] = ref[i] * fft_direct_out[i];
    }
  }
  
  fftw_execute(fft_pinvers);
  
  //Locate the max correlation
  double dMax = 0.0;
  int lag = 0;
  for( int i = 0; i < FFT_Size; i++)
  {
    if(fft_inverse_out[i] > dMax)
    {
      dMax = fft_inverse_out[i];
      lag = i;
    }
  }
  
  if( lag >= FFT_Size/2)
  {
    lag = FFT_Size - lag; 
  }
  else
  {
    lag = -lag;
  }
  
  return lag;
}

NumericVector LabelFreeAlign::getHannWindow()
{
  NumericVector hannWin(WinLength);
  memcpy(hannWin.begin(), HannWindow, sizeof(double)*WinLength);
  return hannWin;
}

NumericVector LabelFreeAlign::getRefLowFFT()
{
  NumericVector refFft(FFT_Size);
  memcpy(refFft.begin(), fft_ref_low, sizeof(double)*FFT_Size);
  return refFft;
}

NumericVector LabelFreeAlign::getRefHighFFT()
{
  NumericVector refFft(FFT_Size);
  memcpy(refFft.begin(), fft_ref_high, sizeof(double)*FFT_Size);
  return refFft;
}

//TESTS//////////////////////////////////////////////////////////////////////////////////////////////
/*
// [[Rcpp::export]]
List TestComputeRefAndHannWin(NumericVector refSpectrum)
{
  double refC[refSpectrum.length()];
  memcpy(refC, refSpectrum.begin(), sizeof(double)*refSpectrum.length());
  LabelFreeAlign alngObj(refC, refSpectrum.length());
  return List::create( Named("HannWin") = alngObj.getHannWindow(), Named("RefLow") = alngObj.getRefLowFFT(), Named("RefHigh") = alngObj.getRefHighFFT());
}

//To run this debug function the ZeroPadding method must be set as public
// [[Rcpp::export]]
NumericVector TestZeroPadding(NumericVector x, bool rev)
{
  const int NewLength = (int)pow(2.0, std::ceil(log2((double)x.length())));
  double refC[NewLength];
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length());
  
  alngObj.ZeroPadding(refC, rev,  NewLength, x.length());
    
  NumericVector y(NewLength);
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y;
}

//To run this debug function the TestTimeWindow method must be set as public
// [[Rcpp::export]]
NumericVector TestTimeWindow(NumericVector x, bool bHigh)
{
  double refC[x.length()];
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length());
  
  alngObj.TimeWindow(refC, bHigh);
  
  NumericVector y(x.length());
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y; 
}

//To run this debug function the FourierLinerScaleShift method must be set as public
// [[Rcpp::export]]
NumericVector TestFourierLinerScaleShift(NumericVector x, double scaling,  double shift)
{
  double refC[x.length()];
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length());
  
  alngObj.FourierLinerScaleShift(refC, scaling, shift);
  
  NumericVector y(x.length());
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y;
}

//To run this debug function the ##### method must be set as public
// [[Rcpp::export]]
double TestFourierBestCor(NumericVector ref, NumericVector x, bool bRefLow)
{
  double refC[ref.length()];
  double xC[x.length()];
  memcpy(refC, ref.begin(), sizeof(double)*ref.length());
  memcpy(xC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, ref.length());
  
  NumericVector ref_ptr;
  if(bRefLow)
  {
    ref_ptr = alngObj.getRefLowFFT();
  }
  else
  {
    ref_ptr = alngObj.getRefHighFFT();
  }
  
  return (double)alngObj.FourierBestCor(xC, ref_ptr.begin());
}

//To run this debug function the ##### method must be set as public
// [[Rcpp::export]]
NumericVector AlignSpectrumToReference(NumericVector ref, NumericVector x)
{
  double refC[ref.length()];
  double xC[x.length()];
  memcpy(refC, ref.begin(), sizeof(double)*ref.length());
  memcpy(xC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, ref.length());
  LabelFreeAlign::TLags lags = alngObj.AlignSpectrum(xC);
  Rcout<<"Lag low = "<<lags.lagLow<<" Lag high = "<<lags.lagHigh<<"\n";
  NumericVector y(x.length());
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y;
}
*/