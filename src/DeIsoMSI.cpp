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

Deisotoper::Deisotoper(IsoDef imgRuninfo, double** PeakMatrix, double* ProMassAxis, double* RawMassAxis)
{
  pPeakOrder = new int[imgRuninfo.massPeaks];
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    pPeakOrder[i] = 0;
  }
  
  pRawPeakOrder = new int[imgRuninfo.massPeaks];
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    pRawPeakOrder[i] = 0;
  }
  
  pTPI = new double[imgRuninfo.massPeaks];
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    pTPI[i] = 0;
  }
  
  PeakFilter(PeakMatrix,imgRuninfo, ProMassAxis, RawMassAxis);
  
  pNumCandidates = new int[imgRuninfo.massPeaks];
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    pNumCandidates[i] = 0;
  }
}

Deisotoper::~Deisotoper()
{
  delete[] pPeakOrder;
  delete[] pRawPeakOrder;
  delete[] pTPI;
  delete[] pNumCandidates;
}

void Deisotoper::PeakFilter(double** PeakMatrix, IsoDef imgRuninfo, double* ProMassAxis, double* RawMassAxis)
{
  //Fill the total pixel intensity vector
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
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
  
  for(int j = 0; j < imgRuninfo.massPeaks; j++)
  {
    tmp = 0;
    indx = 0;
    for(int i = 0; i < imgRuninfo.massPeaks; i++)
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
  
  //Generate the raw peak order vector
  for(int j = 0; j < imgRuninfo.massPeaks; j++)
  {
    tmp = RawMassAxis[imgRuninfo.massChannels - 1];
    for(int k = 0; k < imgRuninfo.massChannels; k++)
    {
      if(std::abs(ProMassAxis[pPeakOrder[j]] - RawMassAxis[k]) > tmp)
      {
        pRawPeakOrder[j] = k - 1;
        break;
      }
      
      if(std::abs(ProMassAxis[pPeakOrder[j]] - RawMassAxis[k]) < tmp)
      {
        tmp = std::abs(ProMassAxis[pPeakOrder[j]] - RawMassAxis[k]);
      }
    }
  }
  
}

NumericMatrix Deisotoper::CandidateFinder(double* ProMassAxis, double* RawMassAxis, IsoDef imgRuninfo)
{
  double CentralIsoMass = 0;
  double LeftIsoMass = 0;
  double RightIsoMass = 0;
  double tmp = 0;
  int cnt = 0;
  int RawMassAxisPointer = 0;
  
  //This loop finds the number of candidates for each peak
  for(int i = 0; i < imgRuninfo.massPeaks; i++)    
  {
    tmp = RawMassAxis[imgRuninfo.massChannels - 1];
    CentralIsoMass = ProMassAxis[pPeakOrder[i]] + (imgRuninfo.CrntNumIso*1.0033548378);
    
    for(int k = pRawPeakOrder[i]; k < imgRuninfo.massChannels; k++)
    {
      if(std::abs(CentralIsoMass - RawMassAxis[k]) < tmp)
      {
        tmp = std::abs(CentralIsoMass - RawMassAxis[k]);
      } else
      {
        RawMassAxisPointer = k - 1;
        break;
      }
    }
    
    LeftIsoMass  = RawMassAxis[RawMassAxisPointer - imgRuninfo.scansT];
    RightIsoMass = RawMassAxis[RawMassAxisPointer + imgRuninfo.scansT];
    
    cnt = 0;
    for(int j = 0; j < imgRuninfo.massPeaks; j++)
    {
      if((pPeakOrder[i] + j) >= imgRuninfo.massPeaks)
      {
        break;
      }
      
      if( ((ProMassAxis[pPeakOrder[i] + j]) > (LeftIsoMass)) && 
          ((ProMassAxis[pPeakOrder[i] + j]) < (RightIsoMass)) )  
      {
        cnt = cnt + 1;
      }
      
      if( (ProMassAxis[pPeakOrder[i] + j]) > (RightIsoMass) )
      {
        break;
      }
    }
    pNumCandidates[i] = cnt;
  }
  
  
  //Building the CanMass matrix
  maxCan = 0;
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    if(maxCan < pNumCandidates[i])
    {
      maxCan = pNumCandidates[i];
    }
  }
  
  //Matrix containing the candidate masses for each selected peak
  NumericMatrix CanMatrix(imgRuninfo.massPeaks, maxCan+1); 
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    for(int j = 1; j <= maxCan; j++)
    {
      CanMatrix(i,j) = 0; 
    }
  }
  
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    CanMatrix(i,0) = pPeakOrder[i]; 
  }
  
  //Computing the candidates peak indexes of each mass
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    tmp = RawMassAxis[imgRuninfo.massChannels - 1];
    CentralIsoMass = ProMassAxis[pPeakOrder[i]] + (imgRuninfo.CrntNumIso*1.0033548378);
    
    for(int k = pRawPeakOrder[i]; k < imgRuninfo.massChannels; k++)
    {
      if(std::abs(CentralIsoMass - RawMassAxis[k]) > tmp)
      {
        RawMassAxisPointer = k - 1;
        break;
      }
      
      if(std::abs(CentralIsoMass - RawMassAxis[k]) < tmp)
      {
        tmp = std::abs(CentralIsoMass - RawMassAxis[k]);
      }
    }
    
    LeftIsoMass  = RawMassAxis[RawMassAxisPointer - imgRuninfo.scansT];
    RightIsoMass = RawMassAxis[RawMassAxisPointer + imgRuninfo.scansT];
    
    cnt = 1; 
    for(int j = 0; j < imgRuninfo.massPeaks; j++)
    {
      if((pPeakOrder[i] + j) >= imgRuninfo.massPeaks)
      {
        break;
      }
      
      if( ((ProMassAxis[pPeakOrder[i] + j]) > (LeftIsoMass)) && 
          ((ProMassAxis[pPeakOrder[i] + j]) < (RightIsoMass)) ) 
      {
        CanMatrix(i,cnt) = (pPeakOrder[i] + j);
        cnt = cnt + 1;
      }
      
      if((ProMassAxis[pPeakOrder[i] + j]) > (RightIsoMass))
      {
        break;
      }
    }
  }
  return CanMatrix;
}

double* Deisotoper::ScoreCalculator(int* test, int NumCan, double** PeakMatrix, IsoDef imgRuninfo, double* result, double* ProMassAxis)
{
  double x[imgRuninfo.numPixels];  //array containing the image of the principal mass 
  double y[imgRuninfo.numPixels];  //array containing the image of the candidate mass
  double A, B, SStot, SSres, m_slope, t_slope, tmpM, tmpI, numC, numCFact, CrntNumIsoFctrl = 0;
  
  for(int i = 0; i < NumCan; i++)
  {
    result[i] = 0;
  }
  
  for(int i = 0; i < imgRuninfo.numPixels ; i++)
  {
    x[i] = PeakMatrix[i][test[0]];  //Monoisotopic peak
  }
  
  for(int i = 0; i < NumCan; i++)
  {
    A = 0;
    B = 0;
    SStot = 0;
    SSres = 0;
    
    for(int j = 0; j < imgRuninfo.numPixels ; j++)
    {
      y[j] = PeakMatrix[j][test[i+1]];  //Isotopic peak
    }
    
    //Model slope
    for(int j = 0; j < imgRuninfo.numPixels; j++)  
    {
      A += (y[j] * x[j]);
      B += (x[j] * x[j]);
    }
    m_slope = A/B;  //Isotope ratio from the data
    
    if (m_slope < 0) //If the slope is negative then is not a possible candidate
    {
      tmpM = 0;
      tmpI = 0;
    }
    else
    { 
      //***********************************//Morphology score (R2)//*******************************************************//
      for(int j = 0; j < imgRuninfo.numPixels ; j++)
      {
        SStot += (y[j]) * (y[j]);  //Total sum of squares: the sum of the squares of the difference of the dependent variable and its mean
      }
      for(int j = 0; j < imgRuninfo.numPixels ; j++)
      {
        SSres += (y[j] - (m_slope*x[j])) * (y[j] - (m_slope*x[j])); //Sum of squared residuals: the sum of the squares of the residuals
      }
      tmpM = 1 - (SSres/SStot);  //Definition of R2 
      
      //***********************************//Intensity Score//*******************************************************//
      
      
      numC  = nearbyint((0.076*ProMassAxis[test[0]] - 12)); //theoretical number of Carbons extracted from the Human metabolome database model
      CrntNumIsoFctrl = 1;
      for(int j = 0; j < imgRuninfo.CrntNumIso; j++)
      {
        CrntNumIsoFctrl *= (j + 1);
      }
      
      numCFact = 1;
      for(int j = 0; j < imgRuninfo.CrntNumIso; j++)
      {
        numCFact *= (numC - j);
      }
      
      t_slope = (numCFact * pow(0.0107, imgRuninfo.CrntNumIso)) / (CrntNumIsoFctrl * pow(0.9893, imgRuninfo.CrntNumIso));     //Computing the theoretical ratio using the isotope distribution formula 
      
      tmpI = (t_slope - m_slope);                // The score is computed using the following formula:
      tmpI = tmpI/0.5;                           //                y = e^(-(x/0.5)^2)
      tmpI = -(tmpI * tmpI);
      tmpI = std::exp(tmpI); 
    }
    
    //Returning the scores
    result[5*i] = tmpM * tmpI;
    result[5*i + 1] = tmpM; 
    result[5*i + 2] = tmpI;
    result[5*i + 3] = m_slope;
    
    if(imgRuninfo.CrntNumIso == 1)
    {
      result[5*i + 4] = m_slope*0.9893/0.0107; //Estimated number of carbons extracted from the first isotope distribution
    }
    
  }
  return result;
}

List Deisotoper::MatrixAnnotator(double **PeakMatrix, IsoDef imgRuninfo, NumericMatrix CanMatrix, double* ProMassAxis)
{
  List AnnMatrix(imgRuninfo.massPeaks);   // List containing the results of the annotation process
  int* CandidateRow;
  CandidateRow = new int[imgRuninfo.MaxCan];
  double* scores;
  scores = new double[imgRuninfo.MaxCan * 5];
  NumericMatrix TMP(6 ,imgRuninfo.MaxCan); 
  
  for(int i = 0; i < 5*imgRuninfo.MaxCan; i++)
  {
    scores[i] = 0;
  }
  
  //For each peak & candidates, the scores are calculated and annotated in the matrix
  for (int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    if(pNumCandidates[i] > 0)
    {
      for(int j = 0; j <= pNumCandidates[i]; j++)
      {
        CandidateRow[j] = CanMatrix(i,j);
      }
      
      scores = ScoreCalculator(CandidateRow, pNumCandidates[i], PeakMatrix, imgRuninfo, scores, ProMassAxis);
      
      for(int k = 0; k < pNumCandidates[i]; k++)
      {
        TMP(0,k) = ProMassAxis[CandidateRow[k + 1]];
        TMP(1,k) = scores[5*k];
        TMP(2,k) = scores[(5*k) + 1];
        TMP(3,k) = scores[(5*k) + 2];
        TMP(4,k) = scores[(5*k) + 3];
        TMP(5,k) = scores[(5*k) + 4];
      }
      
      for(int k = pNumCandidates[i]; k < imgRuninfo.MaxCan; k++)
      {
        TMP(0,k) = -1;
        TMP(1,k) = -1;
        TMP(2,k) = -1;
        TMP(3,k) = -1;
        TMP(4,k) = -1;
        TMP(5,k) = -1;
      }
      
      AnnMatrix[i] = TMP(Range(0,5),Range(0,pNumCandidates[i]-1));
    }
  }
  
  delete[] scores;
  delete[] CandidateRow;
  return AnnMatrix;
}


List Deisotoper::Run(double **PeakMatrix, IsoDef imgRuninfo, double* ProMassAxis ,double* RawMassAxis)
{
  List Result(imgRuninfo.NumIso + 1);
  List PeakResults;
  NumericMatrix CandidatesMatrix;
  NumericVector NUMCAN(imgRuninfo.massPeaks);
  
  for(int i = 0; i < imgRuninfo.NumIso; i++)
  {
    imgRuninfo.CrntNumIso = i + 1;
    CandidatesMatrix = CandidateFinder(ProMassAxis, RawMassAxis, imgRuninfo);
    
    imgRuninfo.MaxCan = maxCan;
    
    PeakResults = MatrixAnnotator(PeakMatrix, imgRuninfo, CandidatesMatrix, ProMassAxis);
    Result[i + 1] = PeakResults;
  }
  
  NumericVector PeakMassesVector(imgRuninfo.massPeaks);
  for(int i = 0; i < imgRuninfo.massPeaks; i++ )
  {
    PeakMassesVector[i] = ProMassAxis[pPeakOrder[i]];
  }
  
  Result[0] = PeakMassesVector;
  
  return Result;
}


// [[Rcpp::export]]
Rcpp::List PeakSelectorC(int massPeaks, int massChannels, int numPixels, int NumIso, NumericMatrix PeakMtx, NumericVector massVec, NumericVector massChanVec, int scansT)
{
  Deisotoper::IsoDef imgRuninfo;
  imgRuninfo.massPeaks = massPeaks;
  imgRuninfo.massChannels = massChannels;
  imgRuninfo.numPixels = numPixels;
  imgRuninfo.NumIso = NumIso;
  imgRuninfo.scansT = scansT;
  
  double** pMatrix;
  pMatrix = new double*[imgRuninfo.numPixels];
  for(int i = 0; i < imgRuninfo.numPixels; i++)
  {
    pMatrix[i] = new double[imgRuninfo.massPeaks];
    for(int j = 0; j < imgRuninfo.massPeaks; j++)
    {
      pMatrix[i][j] = PeakMtx(i,j);
    }
  }
  
  double* ProMassAxis;    //mass axis of the peak matrix
  ProMassAxis = new double[imgRuninfo.massPeaks];
  for(int i = 0; i < imgRuninfo.massPeaks; i++)
  {
    ProMassAxis[i] = massVec[i];
  }
  
  double* RawMassAxis;    //mass axis of the massChannels
  RawMassAxis = new double[imgRuninfo.massChannels];
  for(int i = 0; i < imgRuninfo.massChannels; i++)
  {
    RawMassAxis[i] = massChanVec[i];
  }
  
  //Program call
  Deisotoper myDeisotoper(imgRuninfo, pMatrix, ProMassAxis, RawMassAxis);
  List result = myDeisotoper.Run(pMatrix,imgRuninfo, ProMassAxis, RawMassAxis);
  
  delete[] RawMassAxis;
  
  delete[] ProMassAxis;
  
  for(int i = 0; i < imgRuninfo.numPixels; i++)
  {
    delete[] pMatrix[i];
  }
  delete[] pMatrix;
  
  return result;
}



