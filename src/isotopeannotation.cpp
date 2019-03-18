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
Deisotoper::Deisotoper(IsoDef *myimgRuninfo, NumericMatrix PeakMatrix, NumericVector MatrixAxis, NumericVector ImageAxis)
{
  imgRuninfo = myimgRuninfo;
  
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
  
  massErrorCurve = new double[imgRuninfo->massChannels-1];  //Mass error curve
  for(int i = 0; i < imgRuninfo->massChannels-1; i++)
  {
    massErrorCurve[i] = ((ImageAxis[i+1]-ImageAxis[i])*1000000)/ImageAxis[i+1];
  }
  
  SortPeaksByIntensity(); //Sorts the peaks by intensity
  MatrixToImageAxis(); 
  
  pNumCandidates = new int[imgRuninfo->massPeaks]; //vector containing the number of peak candidates for each selected mass
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    pNumCandidates[i] = 0;
  }
}


//Constructor without Image axis
Deisotoper::Deisotoper(IsoDef *myimgRuninfo, NumericMatrix PeakMatrix, NumericVector MatrixAxis)
{
  imgRuninfo = myimgRuninfo;
  
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
  
  SortPeaksByIntensity(); //Selects the peaks by intensity and points each preocessed mass to the raw mass axis
  
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
  delete[] massErrorCurve;
  
  for(int i; i < imgRuninfo->numPixels; i++)
  {
    delete[] pMatrix[i];
  }
  delete[] pMatrix;
}


//Sorts the peaks by intensity 
void Deisotoper::SortPeaksByIntensity()
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


//Finds the image axis scans closer to the matrix axis masses
void Deisotoper::MatrixToImageAxis()
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
NumericMatrix Deisotoper::CandidateFinder(NumericVector PeaksToCheck)
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


//Gets the tolerance from the mass error curve calculated during the constructor
double Deisotoper::getToleranceFromCurve(int imgIndex)
{
  return(massErrorCurve[imgIndex]);
}


//Computes all the scores over a candidates row 
double* Deisotoper::ScoreCalculator(int* test, int NumCan, double* result, double last_m_slope, int peak_number)
{
  double x[imgRuninfo->numPixels];  //array containing the image of the principal mass 
  double y[imgRuninfo->numPixels];  //array containing the image of the candidate mass
  double A = 0, B = 0, y_mean = 0, x_mean = 0, SStot, SSres, m_slope, intercept;
  double tmpM, tmpI,tmpP, ModCA, CA = 0,ppm = 0 ,maxppm;
  
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
      
      CA = (CrntNumIso == 2) ? (sqrt(1+(8*m_slope*(0.9893/0.0107)*(0.9893/0.0107)))+1)*0.5 : CrntNumIso*((m_slope/last_m_slope)*(0.9893/0.0107) + 1) - 1;

      tmpI  = (ModCA - CA)/(sqrt((1/(log(10)-log(7))))*20);  //Adjusting 20 atoms of difference to score 0.7
      tmpI *= -tmpI;                    
      tmpI  = exp(tmpI); 
    
    //***********************************//Mass ppm Score//******************************************************************//
    
    ppm = (((pMatrixAxis[test[i+1]])-(pMatrixAxis[test[0]]+CrntNumIso*1.0033548378))*1000000)/(pMatrixAxis[test[0]]+CrntNumIso*1.0033548378); //mass error in ppm
    
    maxppm = (!imgRuninfo->ToleranceInScans) ? imgRuninfo->tolerance : getToleranceFromCurve(pImageAxis[pImagePeakOrder[i]]); //maximum valid tolerance
  
    tmpP  =  ppm/(sqrt((1/log(2)))*maxppm); //Adjusting the maximum tolerance to score 0.5
    tmpP *= -tmpP;
    tmpP  = exp(tmpP);
    
    //***********************************//Results//*******************************************************//  
    
      //Returning the scores
      result[(7*i) + 0] = ppm;                      //ppm error
      result[(7*i) + 1] = ((tmpM + tmpI + tmpP)/3); //final score 
      result[(7*i) + 2] = tmpM;                     //morphology score
      result[(7*i) + 3] = tmpI;                     //intensity score
      result[(7*i) + 4] = tmpP;                     //mass error score
      result[(7*i) + 5] = m_slope;                  //model slope  
      result[(7*i) + 6] = CA;                       //number of C atmos

  }

  return result;
}


//Combines the information of the CandidteFinder and the ScoreCalculator
List Deisotoper::MatrixAnnotator(NumericMatrix CanMatrix, NumericVector PeaksToCheck)
{
  List AnnMatrix(imgRuninfo->massPeaks);   // List containing the results of the annotation process
  int* CandidateRow;  //pointer to a vector containing the index of the masses to be tested
  CandidateRow = new int[maxCan];
  double* scores;
  scores = new double[maxCan*7];
  NumericMatrix TMP(10 ,maxCan); 
  
  for(int i = 0; i < 7*maxCan; i++)
  {
    scores[i] = 0;
  }
  
  //For each peak & candidates, the scores are calculated and annotated in the matrix
  for (int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    if((PeaksToCheck[i] > 0) & (pNumCandidates[i] > 0)) //If we need to check the peak and it has candidates , extract the candidates matrix row.
    {
      for(int j = 0; j <= pNumCandidates[i]; j++) 
      {
        CandidateRow[j] = CanMatrix(i,j);
      }
                                
      scores = ScoreCalculator(CandidateRow, pNumCandidates[i], scores, PeaksToCheck[i], i);
      
      for(int k = 0; k < pNumCandidates[i]; k++)
      {
        TMP(0,k) = pMatrixAxis[CandidateRow[k + 1]];  //M+N mass
        TMP(1,k) = scores[(7*k) + 0];                 //ppm error
        TMP(2,k) = CandidateRow[0] + 1;               //M+0 mass index // +1 in order to addapt the indexes to the R format
        TMP(3,k) = CandidateRow[k + 1] + 1;           //M+N mass index
        TMP(4,k) = scores[(7*k) + 1];                 //Final score
        TMP(5,k) = scores[(7*k) + 2];                 //Morphology score
        TMP(6,k) = scores[(7*k) + 3];                 //Intensity score
        TMP(7,k) = scores[(7*k) + 4];                 //Mass error score
        TMP(8,k) = scores[(7*k) + 5];                 //Slope
        TMP(9,k) = scores[(7*k) + 6];                 //Number of C atoms
      }
      
      for(int k = pNumCandidates[i]; k < maxCan; k++)   //This columns will be removed later on
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
        TMP(9,k) = -1; 
      }
      AnnMatrix[i] = TMP(Range(0,9),Range(0,pNumCandidates[i]-1));
    }
  }
  
  delete[] scores;
  delete[] CandidateRow;
  return AnnMatrix;
}


//Program executioning function
List Deisotoper::Run()
{
  List Result(imgRuninfo->numIso + 2);  //List containing the output results
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
  double Hscore = 0; //Highest score at that time
  for(int i = 0; i < imgRuninfo->numIso; i++)  //
  {
    CrntNumIso = i + 1;  //Current searched isotope number
    NumericVector PksTChk = PeaksToCheck(_,i); //Current PeaksToCheck, contains information of the previous stage.

    //Matrix filled with the candidates for each monoisotopic peak
    NumericMatrix CandidatesMatrix = CandidateFinder(PksTChk);  //Generates a matrix that contains all the candidats to be isotopes for each peak mass
    //MaxCan = maxCan;   //Maximum number of candidates that a peak has at the current stage

    PeakResults = MatrixAnnotator(CandidatesMatrix, PksTChk);  //Computes the scores for each candidate in the candidate matrix that are 1st isotopes or have a positive score in the previous stage.
    Result[i + 1] = PeakResults; //Results[0] containts the mass peak masses vector
  
    //TODO CREATE A FUNCTION THAT PACKAGES THIS SHIT!!
    for(int j = 0; j < imgRuninfo->massPeaks; j++) //Filling the PeksToCheck vector with the new results
    {
      if((PksTChk[j] > 0) & (pNumCandidates[j] > 0)) //Check if the tests has been done and if there are any candidates
      {
        NumericMatrix TestResults = PeakResults[j];
        for(int k = 0; k < pNumCandidates[j]; k++)
        {
          if((TestResults(4,k) > imgRuninfo->scoreThreshold) & (TestResults(4,k) > Hscore)) //Checks if the candidate has passed the test
          {
            PeaksToCheck(j,CrntNumIso) = TestResults(8,k);
            Hscore = TestResults(8,k);
            // if(k == 0)
            // {
            //   PeaksToCheck(j,CrntNumIso) = TestResults(8,k);
            // } 
            //   else
            //   {
            //     for(int l = 0; l < k; l++)
            //     {
            //       if((std::abs(TestResults(3,k))) < (std::abs(TestResults(3,l)))) 
            //       {
            //         WriteFlag += 1;
            //       }
            //     }
            //     if(WriteFlag == k)
            //     {
            //       PeaksToCheck(j,CrntNumIso) = TestResults(8,k); 
            //     }
            //     WriteFlag = 0;
            //  }
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
  Deisotoper::IsoDef myimgRuninfo;  
  myimgRuninfo.massPeaks = massPeaks;
  myimgRuninfo.massChannels = massChannels;
  myimgRuninfo.numPixels = numPixels;
  myimgRuninfo.numIso = numIso;
  myimgRuninfo.tolerance = tolerance;
  myimgRuninfo.scoreThreshold = scoreThreshold;
  myimgRuninfo.ToleranceInScans = ToleranceInScans;
  
  List result;

  if(myimgRuninfo.ToleranceInScans)
  {
    Deisotoper myDeisotoper(&myimgRuninfo, PeakMtx, massVec, massChanVec); //Class constructor in scan mode
    result = myDeisotoper.Run();  //Algorithm run
  }else
    {
      Deisotoper myDeisotoper(&myimgRuninfo, PeakMtx, massVec); //Class constructor in ppm mode
      result = myDeisotoper.Run();  //Algorithm run
    }
  
  return result;
}

