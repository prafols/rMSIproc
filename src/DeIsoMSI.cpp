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
  
  PeakFilter(PeakMatrix,imgRuninfo);
  
  pNumCandidates = new int[imgRuninfo.massChannels];
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    pNumCandidates[i] = 0;
  }
}

Deisotyper::~Deisotyper()
{
  delete[] pPeakOrder;
  delete[] pTPI;
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
}

NumericMatrix Deisotyper::CandidateFinder(double* massAxis, IsoDef imgRuninfo)
{
  double MonoIsoMass = 0;
  int cnt = 0;
  
  //This loop finds the number of candidates for each peak
  for(int i = 0; i < imgRuninfo.massChannels; i++)    
  {
    MonoIsoMass = massAxis[pPeakOrder[i]];
    cnt = 0; 
    for(int j = 0; j < imgRuninfo.massChannels; j++)
    {
      if((pPeakOrder[i] + j) >= imgRuninfo.massChannels)
      {
        break;
      }
      
      if( ((massAxis[pPeakOrder[i] + j]) > (MonoIsoMass + 1.003355 - (imgRuninfo.ppmT*MonoIsoMass/10e6))) && 
          ((massAxis[pPeakOrder[i] + j]) < (MonoIsoMass + 1.003355 + (imgRuninfo.ppmT*MonoIsoMass/10e6))) ) 
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
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    if(MaxCan < pNumCandidates[i])
    {
      MaxCan = pNumCandidates[i];
    }
  }
  
  NumericMatrix CanMatrix(imgRuninfo.massChannels,MaxCan+1); //Matrix containing the candidate masses for each selected peak
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    for(int j = 1; j <= MaxCan; j++)
    {
      CanMatrix(i,j) = 0; 
    }
  }
  
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    CanMatrix(i,0) = pPeakOrder[i]; 
  }
  
  //Computing the candidates peak indexes of each mass
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    MonoIsoMass = massAxis[pPeakOrder[i]];
    cnt = 1; 
    for(int j = 0; j < imgRuninfo.massChannels; j++)
    {
      if((pPeakOrder[i] + j) >= imgRuninfo.massChannels)
      {
        break;
      }
      
      if( ((massAxis[pPeakOrder[i] + j]) > (MonoIsoMass + 1.003355 - (imgRuninfo.ppmT*MonoIsoMass/10e6))) && 
          ((massAxis[pPeakOrder[i] + j]) < (MonoIsoMass + 1.003355 + (imgRuninfo.ppmT*MonoIsoMass/10e6))) ) 
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

double* Deisotyper::ScoreCalculator(int* test, int NumCan, double** PeakMatrix, int numPixels, double* result, double* massAxis)
{
  double x[numPixels];  //array containing the image of the principal mass 
  double y[numPixels];  //array containing the image of the candidate mass
  double A = 0;
  double B = 0;
  double SStot = 0;
  double SSres = 0;
  double m_slope = 0;
  
  for(int i = 0; i < NumCan; i++)
  {
  result[i] = 0;
  }

  for(int i = 0; i < numPixels ; i++)
  {
    x[i] = PeakMatrix[i][test[0]];
  }
  
  for(int i = 0; i < NumCan; i++)
  {
    A = 0;
    B = 0;
    SStot = 0;
    SSres = 0;
    
    for(int j = 0; j < numPixels ; j++)
    {
      y[j] = PeakMatrix[j][test[i+1]];
    }
    
    //Model slope
    for(int j = 0; j < numPixels; j++)  
    {
      A += (y[j] * x[j]);
      B += (x[j] * x[j]);
    }
    m_slope = A/B;
    
    if (m_slope < 0)
    {
      result[i] = 0;
    }
      else
      {
      //Morphology score (R2)
      for(int j = 0; j < numPixels ; j++)
      {
        SStot += (y[j]) * (y[j]);  
      }
      for(int j = 0; j < numPixels ; j++)
      {
        SSres += (y[j] - (m_slope*x[j])) * (y[j] - (m_slope*x[j]));  
      }
      result[i] = 1 - (SSres/SStot);  
      
      //multiplying both scores
      result[i] = result[i] * SlopeScore(m_slope, massAxis[test[0]]); //Calls the slope score with the model score in result
      }
    }
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
                                  // TODO: el ratio per el segon isotop
  return score;
}

NumericMatrix Deisotyper::MatrixAnnotator(double **PeakMatrix, IsoDef imgRuninfo, NumericMatrix CanMatrix, double* massAxis)
{
  NumericMatrix AnnMatrix(imgRuninfo.massChannels, MaxCan + MaxCan + 1);
  
  int* CandidateRow;     
  double* result;
  CandidateRow = new int[MaxCan];
  result = new double[MaxCan];
  int massPointer;
  
  for(int i = 0; i < MaxCan; i++)
  {
    result[i] = 0;
  }
  
  //Copying the candidates matrix into the results matrix
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    for(int j = 0; j <= MaxCan; j++)
    {
      massPointer = CanMatrix(i,j);
      if (massPointer > 0)
      {
        AnnMatrix(i,j) = massAxis[massPointer]; //aqui es pot escriure l'index de la columna enlloc de la massa     
      }
    }
  }
  
  //For each peak & candidates, the scores are calculated and annotated in the matrix
  for (int i = 0; i < imgRuninfo.massChannels; i++)
  {
    if(pNumCandidates[i] > 0)
    {
      for(int j = 0; j <= pNumCandidates[i]; j++)
      {
        CandidateRow[j] = CanMatrix(i,j);
      }
      
      result = ScoreCalculator(CandidateRow, pNumCandidates[i], PeakMatrix, imgRuninfo.numPixels, result, massAxis);
      
      for(int j = 0; j < pNumCandidates[i]; j++)
      {
        AnnMatrix(i,  MaxCan + (j + 1) ) = result[j];  
      }
    }
  }

  delete[] result;
  delete[] CandidateRow;
  return AnnMatrix;
}

NumericMatrix Deisotyper::Run(double **PeakMatrix, IsoDef imgRuninfo, double* massAxis)
{
  //return MatrixAnnotator(PeakMatrix, imgRuninfo, CandidateMatrixCheck(PeakFinder(massAxis, imgRuninfo),imgRuninfo), massAxis);
  return MatrixAnnotator(PeakMatrix, imgRuninfo, CandidateFinder(massAxis, imgRuninfo), massAxis);
  
}

// [[Rcpp::export]]
Rcpp::NumericVector PeakSelectorC(int massChannels, int numPixels, int NumIso, NumericMatrix PeakMtx, NumericVector massVec, int ppmT)
{
  Deisotyper::IsoDef imgRuninfo;
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



