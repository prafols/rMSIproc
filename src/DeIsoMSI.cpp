/*************************************************************************
 * 
 *     Copyright (C) 2018 Lluc Sementé Fern?ndez
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

#include "DeIsoMSI.h"
using namespace Rcpp;

Deisotyper::Deisotyper(IsoDef imgRuninfo, double** PeakMatrix)
{
  pPeakOrder = new int[imgRuninfo.massChannels];
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    pPeakOrder[i] = 0;
  }
  
  pTPI = new double[imgRuninfo.massChannels];
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    pTPI[i] = 0;
  }
  
  pAnnotation = new int[imgRuninfo.massChannels];
  
  
  PeakFilter(PeakMatrix,imgRuninfo);
  
  pNumCandidates = new int[FPOlength];
  for(int i = 0; i < FPOlength; i++)
  {
    pNumCandidates[i] = 0;
  }
}


Deisotyper::~Deisotyper()
{
  delete[] pPeakOrder;
  delete[] pTPI;
  delete[] pAnnotation;
  delete[] pNumCandidates;
}
 
 
void Deisotyper::PeakFilter(double** PeakMatrix, IsoDef imgRuninfo)
{
  //Fill the total pixel intensity vector
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    for(int j = 0; j < imgRuninfo.numPixels; j++)
    {
      pTPI[i] += PeakMatrix[j][i];
    }
    pTPI[i] /= imgRuninfo.numPixels;
  }
  
  //Generate the peak order vector
  int indx = 0;    //index of the ordered data 
  double tmp = 0; //temporal maximum value of pTPI  
  int flag = 0;   //index repeat flag
  
  for(int j = 0; j < imgRuninfo.massChannels; j++)
  {
    tmp = 0;
    indx = 0;
    for(int i = 0; i < imgRuninfo.massChannels; i++)
    {
      if(pTPI[i]>tmp)
      {
        flag = 0;
        for(int k = 0; k <= j; k++ )
        {
          if(pPeakOrder[k]==i)
          {
            flag = 1;
          }
        }
        if(flag != 1)
        {
          tmp = pTPI[i];
          indx = i;
        }
      }
    }
    pPeakOrder[j] = indx;
  }
  
  //Selecting the peaks that have enough intensity 
  FPOlength = 1;
  for(int i = 0; i < imgRuninfo.massChannels; i++)
    {
      if (pTPI[pPeakOrder[i]] > (pTPI[pPeakOrder[0]]*imgRuninfo.IntPcent))
      {
        FPOlength++;
      }
        else
        {
          break;
        }
    }
}


NumericMatrix Deisotyper::PeakFinder(double* massAxis, IsoDef imgRuninfo)
{
  double MonoIsoMass = 0;
  int cnt = 0;
  
  for(int i = 0; i < FPOlength; i++)    //This loop finds the number of candidates for each peak
  {
    MonoIsoMass = massAxis[pPeakOrder[i]];
    cnt = 0; 
    for(int j = 0; j < FPOlength; j++)
    {
      if((pPeakOrder[i] + j) >= imgRuninfo.massChannels)
      {
        break;
      }
      
      if( ((massAxis[pPeakOrder[i] + j]) > (MonoIsoMass + 1.003355 - (imgRuninfo.ppmT*MonoIsoMass/10e6))) && 
         ((massAxis[pPeakOrder[i] + j]) < (MonoIsoMass + 1.003355 + (imgRuninfo.ppmT*MonoIsoMass/10e6))) && 
          ((pTPI[pPeakOrder[i] + j]) < (pTPI[pPeakOrder[i]])) )
      {
        cnt = cnt + 1;
      }
      
      if((massAxis[pPeakOrder[i] + j]) > (MonoIsoMass + 1.003355 + (imgRuninfo.ppmT*MonoIsoMass/10e6)))
      {
        break;
      }
    }
    pNumCandidates[i] = cnt;
  }
  
  MaxCan = 0;
  //Building the CanMass matrix
  for(int i = 0; i < FPOlength; i++)
  {
    if(MaxCan < pNumCandidates[i])
    {
      MaxCan = pNumCandidates[i];
    }
  }
  
  NumericMatrix CanMatrix(FPOlength,MaxCan+1); //Matrix containing the candidate masses for each selected peak
  for(int i = 0; i < FPOlength; i++)
  {
    for(int j = 1; j <= MaxCan; j++)
    {
      CanMatrix(i,j) = 0; 
    }
  }
  
  for(int i = 0; i < FPOlength; i++)
  {
    CanMatrix(i,0) = pPeakOrder[i]; 
  }
  
  //Computing the candidates peak indexes of each mass
  for(int i = 0; i < FPOlength; i++)
  {
    MonoIsoMass = massAxis[pPeakOrder[i]];
    cnt = 1; 
    for(int j = 0; j < FPOlength; j++)
    {
      if((pPeakOrder[i] + j) >= imgRuninfo.massChannels)
      {
        break;
      }
      
      if( ((massAxis[pPeakOrder[i] + j]) > (MonoIsoMass + 1.003355 - (imgRuninfo.ppmT*MonoIsoMass/10e6))) && 
          ((massAxis[pPeakOrder[i] + j]) < (MonoIsoMass + 1.003355 + (imgRuninfo.ppmT*MonoIsoMass/10e6))) && 
          ((pTPI[pPeakOrder[i] + j]) < (pTPI[pPeakOrder[i]])) )
      {
        CanMatrix(i,cnt) = (pPeakOrder[i] + j);
        cnt = cnt + 1;
      }
      
      if((massAxis[pPeakOrder[i] + j]) > (MonoIsoMass + 1.003355 + (imgRuninfo.ppmT*MonoIsoMass/10e6)))
      {
        break;
      }
    }
  }
  
  return CanMatrix;
}


NumericMatrix Deisotyper::CandidateMatrixCheck(NumericMatrix CanMatrix)
{
  FPOsubstractor = 0;
  
  for(int i = 0; i < FPOlength; i++)    
  {
    for(int j = 0; j < FPOlength; j++)
    {
      for(int k = 1; k <= pNumCandidates[j]; k++)
      {
        if(CanMatrix(i,0) == CanMatrix(j,k))
        {
          FPOsubstractor++;
          for(int l = 0; l <= pNumCandidates[i]; l++)
          {
            CanMatrix(i,l) = 0;
          }
        }
      }
    }
  }
  return CanMatrix;
}


double* Deisotyper::MorphologyTest(int* test, int NumCan, double** PeakMatrix, int numPixels, double* result)
{
  double x[numPixels];  //array containing the image of the principal mass 
  double y[numPixels];  //array containing the image of the candidate mass 
  
  for(int i = 0; i < numPixels ; i++)
  {
    x[i] = PeakMatrix[i][test[0]];
  }
  
  for(int i = 0; i < NumCan; i++)
  {
    for(int j = 0; j < numPixels ; j++)
    {
      y[j] = PeakMatrix[j][test[i+1]];
    }
    result[i] = gsl_stats_correlation(x, 1, y, 1, numPixels);
  }
  return result;
}


double Deisotyper::IntensityTest(int* test, int Candidate, double** PeakMatrix, int numPixels, double* massAxis)
{
  double x[numPixels];  //array containing the image of the principal mass 
  double y[numPixels];  //array containing the image of the candidate mass 
  double cov = 0;
  double sumsq = 0;
  double result = 0;
  
  for(int i = 0; i < numPixels ; i++)
  {
    x[i] = PeakMatrix[i][test[0]];
  }
  
  for(int j = 0; j < numPixels ; j++)
  {
    y[j] = PeakMatrix[j][test[Candidate + 1]];
  }
  
  gsl_fit_mul(x, 1, y, 1, numPixels, &result, &cov, &sumsq);
  result = SlopeScore(result, massAxis[test[0]]);
  
  return result;
}


double Deisotyper::SlopeScore(double m_slope, double mass)
{
  double score = 0;
  double tmp = 0;
  double t_slope, t_mass; 
  double adduct_mass[3];
  adduct_mass[0] = 1;   // Proton
  adduct_mass[1] = 23;  // Sodium
  adduct_mass[2] = 18;  // Potassium
  
  for(int i = 0; i < 3; i++)
  {
    t_mass = mass - adduct_mass[i];
    t_slope = c0 + c1*t_mass;                 // The model is computed using a lineal regression: y = c0 + c1·x
    tmp = (t_slope - m_slope);                // The score is computed using the following formula:
    tmp = tmp/c;                              //                y = e^(-(x/c)^2)
    tmp = -(tmp * tmp);
    tmp = std::exp(tmp);  
    if (tmp > score)
    {
      score = tmp;
    }
  }

  return score;
}


NumericMatrix Deisotyper::MatrixAnnotator(double **PeakMatrix, IsoDef imgRuninfo, NumericMatrix CanMatrix, double* massAxis)
{
  NumericMatrix AnnMatrix(FPOlength, 2*MaxCan + MaxCan + 1);
  
  int* x;     
  double* mTest;
  double iTest = 0;
  x = new int[MaxCan];
  mTest = new double[MaxCan];
  int massPointer;
  
  //Copying the candidates matrix into the results matrix
  for(int k = 0; k < FPOlength; k++)
  {
    for(int l = 0; l <= MaxCan; l++)
    {
      massPointer = CanMatrix(k,l);
      if (massPointer != 0)
      {
        AnnMatrix(k,l) = massAxis[massPointer];        
      }
    }
  }
  
  //For each peak & candidates, the scores are calculated and annotated in the matrix
  for (int i = 0; i < FPOlength; i++)
  {
    if(pNumCandidates[i] != 0)
    {
      for(int j = 0; j <= pNumCandidates[i]; j++)
      {
        x[j] = CanMatrix(i,j);
      }
      
      mTest = MorphologyTest(x, pNumCandidates[i], PeakMatrix, imgRuninfo.numPixels, mTest);
      
      for(int m = 0; m < pNumCandidates[i]; m++)
      {
        if(mTest[m] >= 0.85)
        {
          iTest = IntensityTest(x, m, PeakMatrix, imgRuninfo.numPixels, massAxis);
          if(iTest > 0)
          {
          AnnMatrix(i,  MaxCan + (2*m + 1) ) = mTest[m];            
          AnnMatrix(i,  MaxCan + (2*m + 2) ) = iTest;
          }
        }
      }
    }
  }

  delete[] mTest;
  delete[] x;
  return AnnMatrix;
}


NumericMatrix Deisotyper::Run(double **PeakMatrix, IsoDef imgRuninfo, double* massAxis)
{
  //return CandidateMatrixCheck(PeakFinder(massAxis, imgRuninfo));
  return MatrixAnnotator(PeakMatrix, imgRuninfo, CandidateMatrixCheck(PeakFinder(massAxis, imgRuninfo)), massAxis);
}


// [[Rcpp::export]]
Rcpp::NumericVector PeakSelectorC(int massChannels, int numPixels, double IntPcent, int NumIso, NumericMatrix PeakMtx, NumericVector massVec, int ppmT)
{
  Deisotyper::IsoDef imgRuninfo;
  imgRuninfo.IntPcent = IntPcent;
  imgRuninfo.massChannels = massChannels;
  imgRuninfo.numPixels = numPixels;
  imgRuninfo.NumIso = NumIso;
  imgRuninfo.ppmT = ppmT;
  
  double** pMatrix;
  pMatrix = new double*[imgRuninfo.numPixels];
  for(int i = 0; i < imgRuninfo.numPixels; i++)
  {
    pMatrix[i] = new double[imgRuninfo.massChannels];
    for(int j = 0; j < imgRuninfo.massChannels; j++)
    {
      pMatrix[i][j] = PeakMtx(i,j);
    }
  }
  
  double* massAxis;
  massAxis = new double[imgRuninfo.massChannels];
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    massAxis[i] = massVec[i];
  }

  //Program call
  Deisotyper myDeisotyper(imgRuninfo, pMatrix);
  NumericMatrix res(myDeisotyper.Run(pMatrix,imgRuninfo, massAxis));
  
  
  delete[] massAxis;
  
  for(int i = 0; i < imgRuninfo.numPixels; i++)
  {
    delete[] pMatrix[i];
  }
  delete[] pMatrix;
  
  return res;
}



