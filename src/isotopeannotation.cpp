/*************************************************************************
 * 
 *     Copyright (C) 2018 Lluc Sementé Fernández
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

#include "isotopeannotation.h"
using namespace Rcpp;

//Constructor with Image axis
Deisotoper::Deisotoper(IsoDef *imgRuninfo, NumericMatrix PeakMatrix, NumericVector MatrixAxis, NumericVector ImageAxis)
{
  pMatrix = new double*[imgRuninfo->numPixels];
  for(int i = 0; i < imgRuninfo->numPixels; i++)
  {
    pMatrix[i] = new double[imgRuninfo->massPeaks];
    for(int j = 0; j < imgRuninfo->massPeaks; j++)
    {
      pMatrix[i][j] = PeakMatrix(i,j);
    }
  }
  
  pMatrixPeakOrder = new int[imgRuninfo->massPeaks]; //vector containing the columns ordered by total pixel intensity
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pMatrixPeakOrder[i] = 0;
  }
  
  pImagePeakOrder = new int[imgRuninfo->massPeaks]; //vector containing the index of the ordered columns by TPI in the raw mass axis
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pImagePeakOrder[i] = 0;
  }
  
  pTPI = new double[imgRuninfo->massPeaks];  //pointer to the total pixel intensity vector
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pTPI[i] = 0;
  }
  
  pMatrixAxis = new double[imgRuninfo->massPeaks];  //Matrix axis
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pMatrixAxis[i] = MatrixAxis[i];
  }
  
  pImageAxis = new double[imgRuninfo->massChannels];  //Image axis
  for(int i = 0; i < imgRuninfo->massChannels; i++)
  {
    pImageAxis[i] = ImageAxis[i];
  }
  
  SortPeaksByIntensity(imgRuninfo); //Sorts the peaks by intensity
  MatrixToImageAxis(imgRuninfo); 
  
  pNumCandidates = new int[imgRuninfo->massPeaks]; //vector containing the number of peak candidates for each selected mass
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pNumCandidates[i] = 0;
  }
}


//Constructor without Image axis
Deisotoper::Deisotoper(IsoDef *imgRuninfo, NumericMatrix PeakMatrix, NumericVector MatrixAxis)
{
  pMatrix = new double*[imgRuninfo->numPixels];
  for(int i = 0; i < imgRuninfo->numPixels; i++)
  {
    pMatrix[i] = new double[imgRuninfo->massPeaks];
    for(int j = 0; j < imgRuninfo->massPeaks; j++)
    {
      pMatrix[i][j] = PeakMatrix(i,j);
    }
  }
  
  pMatrixPeakOrder = new int[imgRuninfo->massPeaks]; //vector containing the columns ordered by total pixel intensity
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pMatrixPeakOrder[i] = 0;
  }
  
  pTPI = new double[imgRuninfo->massPeaks];  //pointer to the total pixel intensity vector
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pTPI[i] = 0;
  }
  
  SortPeaksByIntensity(imgRuninfo); //Selects the peaks by intensity and points each preocessed mass to the raw mass axis
  
  pNumCandidates = new int[imgRuninfo->massPeaks]; //vector containing the number of peak candidates for each selected mass
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pNumCandidates[i] = 0;
  }
}


//Destructor
Deisotoper::~Deisotoper()
{
  delete[] pMatrixPeakOrder;
  delete[] pImagePeakOrder;
  delete[] pTPI;
  delete[] pNumCandidates;
  delete[] pMatrixAxis;
  delete[] pImageAxis;
}


//Selects the peaks by intensity and points each preocessed mass to the raw mass axis
void Deisotoper::SortPeaksByIntensity(IsoDef *imgRuninfo)
{
  //Fill the total pixel intensity vector
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    for(int j = 0; j < imgRuninfo->numPixels; j++)
    {
      pTPI[i] += pMatrix[j][i];
    }
    pTPI[i] /= imgRuninfo->numPixels;
  }
  
  //Generate the peak order vector
  int indx = 0;    //index of the ordered data 
  double tmp = 0; //temporal maximum value of pTPI  
  int flag = 0;   //index repeat flag
  
  for(int j = 0; j < imgRuninfo->massPeaks; j++)
  {
    tmp = 0;
    indx = 0;
    for(int i = 0; i < imgRuninfo->massPeaks; i++)
    {
      if(pTPI[i]>tmp)
      {
        flag = 0;
        for(int k = 0; k <= j; k++ )
        {
          if(pMatrixPeakOrder[k]==i)
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
    pMatrixPeakOrder[j] = indx;
  }
}


//Points one mass axis to the other in order to use the tolerance in scans
void Deisotoper::MatrixToImageAxis(IsoDef *imgRuninfo)
{ 
  int tmp = 0;
  //Generate the raw peak order vector
  for(int j = 0; j < imgRuninfo->massPeaks; j++)
  {
    tmp = pImageAxis[imgRuninfo->massChannels - 1];
    for(int k = 0; k < imgRuninfo->massChannels; k++)
    {
      if(std::abs(pMatrixAxis[pMatrixPeakOrder[j]] - pImageAxis[k]) > tmp)
      {
        pImagePeakOrder[j] = k - 1;
        break;
      }
      
      if(std::abs(pMatrixAxis[pMatrixPeakOrder[j]] - pImageAxis[k]) < tmp)
      {
        tmp = std::abs(pMatrixAxis[pMatrixPeakOrder[j]] - pImageAxis[k]);
      }
    }
  }
}


//Finds the number of isotope candidates for each selected peak and creates the candidate matrix containing the mass and his candidates by rows
NumericMatrix Deisotoper::CandidateFinder(IsoDef *imgRuninfo, NumericVector PeaksToCheck)
{
  double CentralIsoMass = 0;  //Central Isotope mass
  double LeftIsoMass = 0;     //Left Isotope mass
  double RightIsoMass = 0;    //Rigth Isotope mass
  double tmp = 0;
  int cnt = 0;
  int ImageAxisPointer = 0;

  //This loop finds the number of candidates for each peak
  for(int i = 0; i < imgRuninfo->massPeaks; i++)    
  {
    if(PeaksToCheck[i] > 0)
    {
      tmp = pImageAxis[imgRuninfo->massChannels - 1]; //the maximum value the peaks can have 

      CentralIsoMass = pMatrixAxis[pMatrixPeakOrder[i]] + (CrntNumIso*1.0033548378); //central isotopic mass recalculated for each peak
      
      for(int k = pImagePeakOrder[i]; k < imgRuninfo->massChannels; k++) //This loop finds the closest index in the raw mass axis for each cenral candidate mass
      {
        if(std::abs(CentralIsoMass - pImageAxis[k]) < tmp)
        {
          tmp = std::abs(CentralIsoMass - pImageAxis[k]);
        } 
        else
        {
          ImageAxisPointer = k - 1;
          break;
        }
      }
      
      LeftIsoMass  = pImageAxis[ImageAxisPointer - imgRuninfo->tolerance];
      RightIsoMass = pImageAxis[ImageAxisPointer + imgRuninfo->tolerance];
      
      cnt = 0;
      for(int j = 0; j < imgRuninfo->massPeaks; j++)
      {
        if((pMatrixPeakOrder[i] + j) >= imgRuninfo->massPeaks)
        {
          break;
        }
        
        if( ((pMatrixAxis[pMatrixPeakOrder[i] + j]) > (LeftIsoMass)) && 
            ((pMatrixAxis[pMatrixPeakOrder[i] + j]) < (RightIsoMass)) )  
        {
          cnt = cnt + 1;
        }
        
        if( (pMatrixAxis[pMatrixPeakOrder[i] + j]) > (RightIsoMass) )
        {
          break;
        }
      }
      pNumCandidates[i] = cnt;
    }
  }
  
  
  //Building the CanMass matrix
  //Finding the maximum nubmer fo candidates a mass has
  maxCan = 0;
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    if(maxCan < pNumCandidates[i])
    {
      maxCan = pNumCandidates[i];
    }
  }
  
  
  //Matrix containing the candidate masses for each selected peak
  NumericMatrix CanMatrix(imgRuninfo->massPeaks, maxCan+1);      
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    for(int j = 0; j < maxCan+1; j++)
    {
      CanMatrix(i,j) = 0; 
    }
  }
  
  
  //Filling the first column of the candidates matrix with the peak masses index
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    CanMatrix(i,0) = pMatrixPeakOrder[i];
  }
  
  
  //Computing the candidates peak indexes of each mass
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    if((PeaksToCheck[i] > 0) & (pNumCandidates[i] > 0))
      {
      tmp = pImageAxis[imgRuninfo->massChannels - 1];
      CentralIsoMass = pMatrixAxis[pMatrixPeakOrder[i]] + (CrntNumIso*1.0033548378);
      
      for(int k = pImagePeakOrder[i]; k < imgRuninfo->massChannels; k++)
      {
        if(std::abs(CentralIsoMass - pImageAxis[k]) > tmp)
        {
          ImageAxisPointer = k - 1;
          break;
        }
        
        if(std::abs(CentralIsoMass - pImageAxis[k]) < tmp)
        {
          tmp = std::abs(CentralIsoMass - pImageAxis[k]);
        }
      }
      
      LeftIsoMass  = pImageAxis[ImageAxisPointer - imgRuninfo->tolerance];
      RightIsoMass = pImageAxis[ImageAxisPointer + imgRuninfo->tolerance];
      
      cnt = 1; 
      for(int j = 0; j < imgRuninfo->massPeaks; j++)
      {
        if((pMatrixPeakOrder[i] + j) >= imgRuninfo->massPeaks)
        {
          break;
        }
        
        if( ((pMatrixAxis[pMatrixPeakOrder[i] + j]) > (LeftIsoMass)) && 
            ((pMatrixAxis[pMatrixPeakOrder[i] + j]) < (RightIsoMass)) ) 
        {
          CanMatrix(i, cnt) = (pMatrixPeakOrder[i] + j);
          cnt++;
        }
        
        if((pMatrixAxis[pMatrixPeakOrder[i] + j]) > (RightIsoMass))
        {
          break;
        }
      }
    }
  }
  return CanMatrix;
}


//Computes the morphology & intensity scores test
double* Deisotoper::ScoreCalculator(int* test, int NumCan,IsoDef *imgRuninfo, double* result, double last_m_slope)
{
  double x[imgRuninfo->numPixels];  //array containing the image of the principal mass 
  double y[imgRuninfo->numPixels];  //array containing the image of the candidate mass
  double A = 0, B = 0, y_mean = 0, x_mean = 0, SStot, SSres, m_slope, intercept, tmpM, tmpI,tmpP, ModCA, CA = 0,ppm = 0;
  
  for(int i = 0; i < NumCan; i++)
  {
    result[i] = 0;
  }
  
  for(int i = 0; i < imgRuninfo->numPixels ; i++)
  {
    x[i] = pMatrix[i][test[0]];  //Monoisotopic peak
  }
  
  for(int j = 0; j < imgRuninfo->numPixels ; j++)
  {
    x_mean += x[j];  //Candidate mass intensity mean
  }
  x_mean = (x_mean/imgRuninfo->numPixels);
  
  //***********************************//Run the test for each candidate//***************************************//
  for(int i = 0; i < NumCan; i++)   
  {
    A = 0;
    B = 0;
    SStot = 0;
    SSres = 0;
    y_mean = 0;
    intercept = 0;
    for(int j = 0; j < imgRuninfo->numPixels ; j++)
    {
      y[j] = pMatrix[j][test[i+1]];  //M+N peak
    }
    
    for(int j = 0; j < imgRuninfo->numPixels ; j++)
    {
      y_mean += y[j];  //Candidate mass intensity mean
    }
    y_mean = (y_mean/imgRuninfo->numPixels);
    
    //***********************************//Lineal model//*******************************************************//
      
      for(int j = 0; j < imgRuninfo->numPixels; j++)  
      {
        A += (x[j] - x_mean)*(y[j] - y_mean);
        B += (x[j] - x_mean)*(x[j] - x_mean);
      }
      m_slope = A/B;  //Isotope ratio from the data
      intercept = y_mean - m_slope*x_mean;
      
    //***********************************//Morphology score (adj. R2)//*******************************************************//
    
      //Total sum of squares: the sum of the squares of the difference of the dependent variable and its mean
      for(int j = 0; j < imgRuninfo->numPixels ; j++)
      {
        SStot += (y[j] - y_mean) * (y[j] - y_mean);     
      }
      
      //Sum of squared residuals: the sum of the squares of the residuals
      for(int j = 0; j < imgRuninfo->numPixels ; j++)
      {
        SSres += (y[j] - (intercept + (m_slope*x[j]))) * (y[j] - (intercept + (m_slope*x[j])));     //TODO probar treure el model i ficar directament el valor de x
      }
      
      //Definition of R2 
      tmpM = 1 - (SSres/SStot); 
      tmpM = tmpM;
      
    //***********************************//Intensity Score//******************************************************************//
    
      //Theoretical number of Carbons extracted from the Human metabolome database model  
      ModCA  = 0.076*pMatrixAxis[test[0]] - 12;
      
      //Number of carbons from the experimental intensity rate
      if(CrntNumIso == 2)
      {
        CA = (std::sqrt(1+(8*m_slope*(0.9893/0.0107)*(0.9893/0.0107)))+1)*0.5;
      }else
        {
          CA = CrntNumIso*((m_slope/last_m_slope)*(0.9893/0.0107) + 1) - 1;   
        }
  
      tmpI = ModCA - CA;                // The score is computed using the following formula:
      tmpI = tmpI/37;                     //                y = e^(-(x/37)^2)
      tmpI = -(tmpI * tmpI);              //TODO adjust intensity score
      tmpI = std::exp(tmpI); 
    
    //***********************************//Mass ppm Score//******************************************************************//
    
    ppm = (((pMatrixAxis[test[i+1]])-(pMatrixAxis[test[0]]+CrntNumIso*1.0033548378))*1000000)/(pMatrixAxis[test[0]]+CrntNumIso*1.0033548378);

    //***********************************//Results//*******************************************************//  
    
      //Returning the scores
      result[(6*i) + 0] = ((tmpM + tmpI)*0.5); //Final score 
      result[(6*i) + 1] = tmpM;                //Morphology score
      result[(6*i) + 2] = tmpI;                //Intensity score
      result[(6*i) + 3] = m_slope;             //Model slope  
      result[(6*i) + 4] = CA;                  //Number of C atoms
      result[(6*i) + 5] = ppm;                 //ppm error
  }

  return result;
}


//Creates the results list for each isotope index
List Deisotoper::MatrixAnnotator(IsoDef *imgRuninfo, NumericMatrix CanMatrix, NumericVector PeaksToCheck)
{
  List AnnMatrix(imgRuninfo->massPeaks);   // List containing the results of the annotation process
  int* CandidateRow;  //pointer to a vector containing the index of the masses to be tested
  CandidateRow = new int[MaxCan];
  double* scores;
  scores = new double[MaxCan*6];
  NumericMatrix TMP(9 ,MaxCan); 
  
  for(int i = 0; i < 6*MaxCan; i++)
  {
    scores[i] = 0;
  }
  
  //For each peak & candidates, the scores are calculated and annotated in the matrix
  for (int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    if((PeaksToCheck[i] > 0) & (pNumCandidates[i] > 0))
    {
      for(int j = 0; j <= pNumCandidates[i]; j++)
      {
        CandidateRow[j] = CanMatrix(i,j);
      }
                                
      scores = ScoreCalculator(CandidateRow, pNumCandidates[i], imgRuninfo, scores, PeaksToCheck[i]);
      
      for(int k = 0; k < pNumCandidates[i]; k++)
      {
        TMP(0,k) = pMatrixAxis[CandidateRow[k + 1]];  //M+N mass
        TMP(1,k) = CandidateRow[0] + 1;               //M+0 mass index // +1 in order to addapt the indexes to the R format
        TMP(2,k) = CandidateRow[k + 1] + 1;           //M+N mass index
        TMP(3,k) = scores[(6*k) + 5];                 //ppm error
        TMP(4,k) = scores[(6*k) + 0];                 //Final score
        TMP(5,k) = scores[(6*k) + 1];                 //Morphology score
        TMP(6,k) = scores[(6*k) + 2];                 //Intensity score
        TMP(7,k) = scores[(6*k) + 3];                 //Slope
        TMP(8,k) = scores[(6*k) + 4];                 //Number of C atoms
      }
      
      for(int k = pNumCandidates[i]; k < MaxCan; k++)
      {
        TMP(0,k) = -1;
        TMP(1,k) = -1;
        TMP(2,k) = -1;
        TMP(3,k) = -1;
        TMP(4,k) = -1;
        TMP(5,k) = -1;
        TMP(6,k) = -1;
        TMP(7,k) = -1;
        TMP(8,k) = -1;
      }
      
      AnnMatrix[i] = TMP(Range(0,8),Range(0,pNumCandidates[i]-1));
    }
  }
  
  delete[] scores;
  delete[] CandidateRow;
  return AnnMatrix;
}


//Program executioning function
List Deisotoper::Run(IsoDef *imgRuninfo)
{
  List Result(imgRuninfo->numIso + 1);  //List containing the output results
  List PeakResults; //Results for each isotope index
  NumericMatrix PeaksToCheck(imgRuninfo->massPeaks, imgRuninfo->numIso); //Matrix that indicates if a ion must undergo the test
                                                                         //A 0 means no check
                                                                         //This matrix is overwritten with the slope of the previous stage if the score is higher than threshold

  //***********************************//Fill the PeaksToCheck matrix//***********************************************// 
  for(int j = 0; j < imgRuninfo->numIso; j++)  
  {
    if(j == 0)
      {
        for(int i = 0; i < imgRuninfo->massPeaks; i++)
        {
          PeaksToCheck(i,j) = 1;  //For the first isotope, all peaks are gonna be checked
        }
      }
        else
        {
          for(int i = 0; i < imgRuninfo->massPeaks; i++)
          {
            PeaksToCheck(i,j) = 0; //By default 0 means no need to check
          } 
        }
  }
 
 
 //***********************************//Algorithm Run (Main structure)//***********************************************//  
  int WriteFlag = 0;
  for(int i = 0; i < imgRuninfo->numIso; i++)  //
  {
    CrntNumIso = i + 1;  //Current searched isotope number
    NumericVector PksTChk = PeaksToCheck(_,i); //Current PeaksToCheck, contains information of the previous stage.

    //Matrix filled with the candidates for each monoisotopic peak
    NumericMatrix CandidatesMatrix = CandidateFinder(imgRuninfo, PksTChk);  //Generates a matrix that contains all the candidats to be isotopes for each peak mass
    MaxCan = maxCan;   //Maximum number of candidates that a peak has at the current stage

    PeakResults = MatrixAnnotator(imgRuninfo, CandidatesMatrix, PksTChk);  //Computes the scores for each candidate in the candidate matrix that are 1st isotopes or have a positive score in the previous stage.
    Result[i + 1] = PeakResults; //Results[0] containts the mass peak masses vector
  
    //TODO CREATE A FUNCTION THAT PACKAGES THIS SHIT!!
    for(int j = 0; j < imgRuninfo->massPeaks; j++) //Filling the PeksToCheck vector with the new results
    {
      if((PksTChk[j] > 0) & (pNumCandidates[j] > 0)) //Check if the tests has been done and if there are any candidates
      {
        NumericMatrix TestResults = PeakResults[j];
        for(int k = 0; k < pNumCandidates[j]; k++)
        {
          if(TestResults(4,k) > imgRuninfo->scoreThreshold) //Checks if the candidate has passed the test
          {
            if(k == 0)
            {
              PeaksToCheck(j,CrntNumIso) = TestResults(7,k);
            } 
              else
              {
                for(int l = 0; l < k; l++)
                {
                  if((std::abs(TestResults(3,k))) < (std::abs(TestResults(3,l)))) 
                  {
                    WriteFlag += 1;
                  }
                }
                if(WriteFlag == k)
                {
                  PeaksToCheck(j,CrntNumIso) = TestResults(7,k); 
                }
                WriteFlag = 0;
              }
            }
          }
        }
      }
    }
  
  
  //***********************************//Adding the masses of the peaks to format the output//***********************************************// 
  NumericVector PeakMassesVector(imgRuninfo->massPeaks); //Vector containing the masses names
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
     PeakMassesVector[i] = pMatrixAxis[pMatrixPeakOrder[i]];
  }
  Result[0] = PeakMassesVector;
  
  
  return Result;
}


// [[Rcpp::export]]
Rcpp::List IsotopeAnnotator(int massPeaks, int massChannels, int numPixels, int numIso,
                            NumericMatrix PeakMtx, NumericVector massVec, NumericVector massChanVec,
                            int tolerance, double scoreThreshold, bool ToleranceInScans)
{
  //Fill the class data structure with the information of the experiment
  Deisotoper::IsoDef imgRuninfo;  
  imgRuninfo.massPeaks = massPeaks;
  imgRuninfo.massChannels = massChannels;
  imgRuninfo.numPixels = numPixels;
  imgRuninfo.numIso = numIso;
  imgRuninfo.tolerance = tolerance;
  imgRuninfo.scoreThreshold = scoreThreshold;
  imgRuninfo.ToleranceInScans = ToleranceInScans;
  List result;
  
  if(imgRuninfo.ToleranceInScans)
  {
    Deisotoper myDeisotoper(&imgRuninfo, PeakMtx, massVec, massChanVec); //Class constructor in scan mode
    result = myDeisotoper.Run(&imgRuninfo);  //Algorithm run
  }else
    {
      Deisotoper myDeisotoper(&imgRuninfo, PeakMtx, massVec); //Class constructor in ppm mode
      result = myDeisotoper.Run(&imgRuninfo);  //Algorithm run
    }
  
  return result;
}

