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

LabelFreeAlign::LabelFreeAlign(double *ref_spectrum, int numOfPoints, bool bilinear, boost::mutex *sharedMutex, int iterations, 
                               double lagRefLow, double lagRefMid, double lagRefHigh,
                               double lagLimitppm, int fftOverSampling, double winSizeRelative ):
dataLength(numOfPoints),
bBilinear(bilinear),
AlignIterations(iterations),
RefLagLow(lagRefLow),
RefLagMid(lagRefMid),
RefLagHigh(lagRefHigh),
FFTScaleShiftOverSampling(fftOverSampling)

{
  WinLength = (int)round(winSizeRelative * (double)dataLength);
  FFT_Size_SH = (int)pow(2.0, std::ceil(log2(dataLength)));

  //Pre-compute the Hanning Windows
  HannWindow = new double[WinLength];
  HannWindowCenter = new double[WinLength];
  for(int i = 0; i < WinLength; i++)
  {
    HannWindow[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i + 1.0) / ( 2.0 * WinLength ) ));
    HannWindowCenter[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i + 1.0) / ( WinLength ) ));
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
  fft_ref_center = new  double[FFT_Size];
  fft_ref_high = new  double[FFT_Size];
  
  ComputeRef(ref_spectrum, BOTTOM_SPECTRUM);
  ComputeRef(ref_spectrum, CENTER_SPECTRUM);
  ComputeRef(ref_spectrum, TOP_SPECTRUM);

  //Lag limits relative to number of points
  lagMax = (lagLimitppm/1.0e6)*(double)dataLength;
  
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
  delete[] fft_ref_center;
  delete[] fft_ref_high;
  delete[] HannWindow;
  delete[] HannWindowCenter;
}

void LabelFreeAlign::ComputeRef(double *data_ref, int spectrumPart)
{
  double *data_ptr;
  switch(spectrumPart)
  {
    case BOTTOM_SPECTRUM:
      data_ptr = fft_ref_low; 
      break;
      
    case CENTER_SPECTRUM:
      data_ptr = fft_ref_center; 
      break;
      
    case TOP_SPECTRUM:
      data_ptr = fft_ref_high;
      break;
  }
  CopyData2Window(data_ref, fft_direct_in, spectrumPart);
  TimeWindow(fft_direct_in,  spectrumPart);
  ZeroPadding(fft_direct_in, spectrumPart == TOP_SPECTRUM, FFT_Size, WinLength);
  fftw_execute(fft_pdirect);
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

void LabelFreeAlign::CopyData2Window(double *data_int, double *data_out, int spectrumPart)
{
  int offset_in;
  int offset_out;
  switch(spectrumPart)
  {
    case BOTTOM_SPECTRUM:
      offset_in = (int)round(RefLagLow * (double)dataLength);
      offset_out = 0;
      break;
      
    case CENTER_SPECTRUM:
      offset_in = (int)round(RefLagMid * (double)dataLength) - WinLength/2;
      offset_out = 0;
      break;
      
    case TOP_SPECTRUM:
      offset_in = (int)round(RefLagHigh * (double)dataLength) - WinLength;
      offset_out = FFT_Size - WinLength;
      break;
  }
  
  //Move pointer offsets to avoid a segfault accesing outside vector
  if( (offset_in + WinLength) >  dataLength)
  {
    offset_in -= (offset_in + WinLength) - dataLength;
  }
  if( offset_in < 0 )
  {
    offset_in = 0;
  }
  
  memcpy(data_out + offset_out, data_int + offset_in, sizeof(double)*WinLength);
}

void LabelFreeAlign::TimeWindow(double *data,  int spectrumPart)
{
  for( int i = 0; i < WinLength; i++)
  {
    
    switch(spectrumPart)
    {
    case BOTTOM_SPECTRUM:
      data[i] *= HannWindow[WinLength - i - 1]; 
      break;
      
    case CENTER_SPECTRUM:
      data[i] *= HannWindowCenter[i];
      break;
      
    case TOP_SPECTRUM:
      data[i + (FFT_Size - WinLength)] *= HannWindow[i]; 
      break;
    }
  }
}

LabelFreeAlign::TLags LabelFreeAlign::AlignSpectrum(double *data )
{
  //Hanning Windowing
  double *topWin_data = new double[FFT_Size];
  double *midWin_data = new double[FFT_Size];
  double *botWin_data = new double[FFT_Size];
  TLags firstLag;
  
  for(int i = 0; i < AlignIterations; i++)
  {
    CopyData2Window(data, topWin_data, TOP_SPECTRUM);
    CopyData2Window(data, midWin_data, CENTER_SPECTRUM);
    CopyData2Window(data, botWin_data, BOTTOM_SPECTRUM);
    
    TimeWindow(topWin_data, TOP_SPECTRUM);
    TimeWindow(midWin_data, CENTER_SPECTRUM);
    TimeWindow(botWin_data, BOTTOM_SPECTRUM);
    
    //Zero-padding 2 improve fft performance
    ZeroPadding(topWin_data, true, FFT_Size, WinLength);
    ZeroPadding(botWin_data, false, FFT_Size, WinLength);
    ZeroPadding(midWin_data, false, FFT_Size, WinLength);
    
    //Get lags
    TLags lags;
    lags.lagLow = (double)FourierBestCor(botWin_data, fft_ref_low);
    lags.lagMid = (double)FourierBestCor(midWin_data, fft_ref_center);
    lags.lagHigh = (double)FourierBestCor(topWin_data, fft_ref_high);
    
    //Limit lags
    lags.lagLow = abs(lags.lagLow) > lagMax ? 0.0 : lags.lagLow;
    lags.lagMid = abs(lags.lagMid) > lagMax ? 0.0 : lags.lagMid;
    lags.lagHigh = abs(lags.lagHigh) > lagMax ? 0.0 : lags.lagHigh; 
    if(i == 0)
    {
      firstLag = lags;
    }
  
    //Aligment constansts mz'(n) = K*mz(n) + Sh
    const double Rl = RefLagLow * (double)dataLength;
    const double Rm = RefLagMid * (double)dataLength;
    const double Rh = RefLagHigh * (double)dataLength;
    
    //Spectra warping constants
    double K1, K2, Sh1, Sh2;
    
    if(bBilinear)
    {
      K1 = (Rm + lags.lagMid - Rl - lags.lagLow)/(Rm - Rl); //New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
      Sh1 = (Rm* lags.lagLow - Rl*lags.lagMid)/(Rm - Rl); //If scaling is performed before shift. New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
      
      K2 = (Rh + lags.lagHigh - Rm - lags.lagMid)/(Rh - Rm); //New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
      Sh2 = (Rh* lags.lagMid - Rm*lags.lagHigh)/(Rh - Rm); //If scaling is performed before shift. New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
      
      //Apply the scaling and shift to original data, after that the pointer to data contains the aligned sptectrum, so copy data before
      double *data2 = new double[dataLength];
      memcpy(data2, data, sizeof(double)*dataLength);
      
      FourierLinerScaleShift(data, K1, Sh1); //Now data will contain the aligned left part of spectrum
      FourierLinerScaleShift(data2, K2, Sh2);
      
      //Merge the two vectors in a single one using the Rm as a center...
      const int iVectorCenterOffset =  dataLength/2; //0.5 es the central reference
      memcpy(data+iVectorCenterOffset, data2+iVectorCenterOffset, sizeof(double)*(dataLength - iVectorCenterOffset));
      delete[] data2;
    }
    else
    {
      K1 = (Rh + lags.lagHigh - Rl - lags.lagLow)/(Rh - Rl); //New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
      Sh1 = (Rh* lags.lagLow - Rl*lags.lagHigh)/(Rh - Rl); //If scaling is performed before shift. New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
      
      //Apply the scaling and shift to original data, after that the pointer to data contains the aligned sptectrum
      FourierLinerScaleShift(data, K1, Sh1);
    }
  }
  
  //Apply a 3 samples moving average to compensate for FFT artifacts
  MovingAverage3Samples(data);
  
  delete[] topWin_data;
  delete[] botWin_data;
  delete[] midWin_data;
    
  return firstLag;
}

void LabelFreeAlign::FourierLinerScaleShift(double *data, double scaling, double shift)
{
  const int NewFFT_Size = FFT_Size_SH + (int)round( (scaling - 1.0) * (double)FFT_Size_SH);
  const int NewFFT_SizeOversampled = NewFFT_Size * FFTScaleShiftOverSampling;
  fftwMtx->lock();
  fftw_complex *fft_odd_in = fftw_alloc_complex(NewFFT_SizeOversampled);
  fftw_complex *fft_odd_out = fftw_alloc_complex(NewFFT_SizeOversampled);
  fftw_plan fft_pOddInvers = fftw_plan_dft_1d(NewFFT_SizeOversampled, fft_odd_in, fft_odd_out, FFTW_BACKWARD, FFTW_ESTIMATE);
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
    
    maxInit = fft_direct_shiftScale_in[i][0] > maxInit ? fft_direct_shiftScale_in[i][0] : maxInit;
  }
  fftw_execute(fft_pdshiftScale);
  
  //Copy data to fft_odd_in and apply scaling and shift
  double arg, Re, Im, auxRe, auxIm;
  for( int i = 0; i < NewFFT_SizeOversampled; i++ )
  {
    //Scaling in FFT domain
    if( i < FFT_Size_SH/2 && i < NewFFT_SizeOversampled/2)
    {
      fft_odd_in[i][0] = fft_direct_shiftScale_out[i][0]; //Copy real part
      fft_odd_in[i][1] = fft_direct_shiftScale_out[i][1]; //Copy imaginary part
    }
    else if( i >= (NewFFT_SizeOversampled - FFT_Size_SH/2))
    {
      fft_odd_in[i][0] = fft_direct_shiftScale_out[i - NewFFT_SizeOversampled + FFT_Size_SH][0]; //Copy real part
      fft_odd_in[i][1] = fft_direct_shiftScale_out[i - NewFFT_SizeOversampled + FFT_Size_SH][1]; //Copy imaginary part
    }
    else
    {
      //Zero padding
      fft_odd_in[i][0] = 0.0;
      fft_odd_in[i][1] = 0.0;
    }
    
    //Shift in FFT domain
    arg = -2.0*M_PI*(1.0/((double)NewFFT_SizeOversampled))*((double)i + 1.0)*round( ((double)FFTScaleShiftOverSampling) * shift);
    Re = cos(arg);
    Im = sin(arg);
    auxRe = fft_odd_in[i][0] * Re - fft_odd_in[i][1]*Im;
    auxIm = fft_odd_in[i][0] * Im + fft_odd_in[i][1]*Re;
    fft_odd_in[i][0] = auxRe;
    fft_odd_in[i][1] = auxIm;
  }
  fftw_execute(fft_pOddInvers);

  //Copy data to output vector downsampling it by averaging values
  double maxEnd = 0.0;
  int lastiUp = -FFTScaleShiftOverSampling/2;
  int iup = 0;
  for( int i = 0; i < dataLength; i++)
  {
    data[i] = 0.0;
    while( iup <  ( lastiUp + FFTScaleShiftOverSampling ))
    {
      data[i] += sqrt(fft_odd_out[iup][0]*fft_odd_out[iup][0] + fft_odd_out[iup][1]*fft_odd_out[iup][1]); //Module
      iup++;
    }
    lastiUp = iup;
    
    //It is not necessary to divide data[i] by FFTScaleShiftOverSampling since I'll compensate for the max.
    maxEnd = data[i] > maxEnd ? data[i] : maxEnd;
  }
  
  //Compensate gain
  if( maxEnd > 0.0 && maxInit > 0.0 )
  {
    for(int i = 0; i < dataLength; i++)
    {
      data[i] /= maxEnd/maxInit;
    }
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

void LabelFreeAlign::MovingAverage3Samples(double *data)
{
  double curr;
  double ant = data[0];
  for( int i = 1; i < (dataLength - 1); i++)
  {
    curr = data[i]; //Current sample backup
    data[i] += data[i+1] + ant;
    data[i] /= 3.0;
    ant = curr;
  }
}

NumericVector LabelFreeAlign::getHannWindow()
{
  NumericVector hannWin(WinLength);
  memcpy(hannWin.begin(), HannWindow, sizeof(double)*WinLength);
  return hannWin;
}

NumericVector LabelFreeAlign::getHannWindowCenter()
{
  NumericVector hannWin(WinLength);
  memcpy(hannWin.begin(), HannWindowCenter, sizeof(double)*WinLength);
  return hannWin;
}

NumericVector LabelFreeAlign::getRefLowFFT()
{
  NumericVector refFft(FFT_Size);
  memcpy(refFft.begin(), fft_ref_low, sizeof(double)*FFT_Size);
  return refFft;
}

NumericVector LabelFreeAlign::getRefCenterFFT()
{
  NumericVector refFft(FFT_Size);
  memcpy(refFft.begin(), fft_ref_center, sizeof(double)*FFT_Size);
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
  boost::mutex mtx;
  double refC[refSpectrum.length()];
  memcpy(refC, refSpectrum.begin(), sizeof(double)*refSpectrum.length());
  LabelFreeAlign alngObj(refC, refSpectrum.length(), false, &mtx);
  return List::create( Named("HannWin") = alngObj.getHannWindow(), Named("HannWinCenter") =  alngObj.getHannWindowCenter(), Named("RefLow") = alngObj.getRefLowFFT(), Named("RefCenter") = alngObj.getRefCenterFFT(),  Named("RefHigh") = alngObj.getRefHighFFT());
}
*/

/*
//To run this debug function the ZeroPadding method must be set as public
// [[Rcpp::export]]
NumericVector TestZeroPadding(NumericVector x, bool rev)
{
  boost::mutex mtx;
  const int NewLength = (int)pow(2.0, std::ceil(log2((double)x.length())));
  double refC[NewLength];
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length(), false, &mtx);
  
  alngObj.ZeroPadding(refC, rev,  NewLength, x.length());
    
  NumericVector y(NewLength);
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y;
}
*/

/*
//To run this debug function the TestTimeWindow method must be set as public
// [[Rcpp::export]]
NumericVector TestTimeWindow(NumericVector x, bool bHigh)
{
  boost::mutex mtx;
  double refC[x.length()];
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length(), false, &mtx);
  
  alngObj.TimeWindow(refC, bHigh);
  
  NumericVector y(x.length());
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y; 
}
*/

/*
//To run this debug function the FourierLinerScaleShift method must be set as public
// [[Rcpp::export]]
NumericVector TestFourierLinerScaleShift(NumericVector x, double scaling,  double shift)
{
  double refC[x.length()];
  boost::mutex mtx;
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length(), false, &mtx);
  
  alngObj.FourierLinerScaleShift(refC, scaling, shift);
  
  NumericVector y(x.length());
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y;
}
*/

/*
//To run this debug function the ##### method must be set as public
// [[Rcpp::export]]
double TestFourierBestCor(NumericVector ref, NumericVector x, bool bRefLow)
{
  boost::mutex mtx;
  double refC[ref.length()];
  double xC[x.length()];
  memcpy(refC, ref.begin(), sizeof(double)*ref.length());
  memcpy(xC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, ref.length(), false, &mtx);
  
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
*/

//To run for debug purposes
// [[Rcpp::export]]
NumericVector AlignSpectrumToReference(NumericVector ref, NumericVector x, bool bilinear = false, 
                                       double lagRefLow = 0.1, double lagRefMid = 0.5, double lagRefHigh = 0.9,
                                       int iterations = 1, double lagLimitppm = 200, int fftOverSampling = 2 )
{
  double *refC = new double[ref.length()];
  double *xC = new double[x.length()];
  
  boost::mutex mtx;
  memcpy(refC, ref.begin(), sizeof(double)*ref.length());
  memcpy(xC, x.begin(), sizeof(double)*x.length());
 
  LabelFreeAlign alngObj(refC, ref.length(), bilinear, &mtx, iterations, lagRefLow, lagRefMid, lagRefHigh, lagLimitppm, fftOverSampling);
  LabelFreeAlign::TLags lags = alngObj.AlignSpectrum(xC);

  Rcout<<"Lag low = "<<lags.lagLow<<" Lag center = "<<lags.lagMid<<" Lag high = "<<lags.lagHigh<<"\n";
  NumericVector y(x.length());
  memcpy(y.begin(), xC, sizeof(double)*x.length());

  delete[] refC;
  delete[] xC;
  return y;
}
