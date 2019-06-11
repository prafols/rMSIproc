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

  if(imgRuninfo->ToleranceInScans)  
  {
    pImagePeakOrder = new int[imgRuninfo->massPeaks]; //vector containing the index of the ordered columns by TPI in the raw mass axis
    for(int i = 0; i < imgRuninfo->massPeaks; i++)
    {
      pImagePeakOrder[i] = 0;
    }
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

  SortPeaksByIntensity(); //Sorts the peaks by intensity

  if(imgRuninfo->ToleranceInScans)  
  {
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
    
    MatrixToImageAxis(); 
  }

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
  delete[] pTPI;
  delete[] pNumCandidates;
  delete[] pMatrixAxis;
  
  if(imgRuninfo->ToleranceInScans)
  {
    delete[] pImagePeakOrder;
    delete[] pImageAxis;
    delete[] massErrorCurve;
  }
  
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

//(Scan mode)Finds the number of isotope candidates for each selected peak and creates the candidate matrix containing the mass and his candidates by rows
NumericMatrix Deisotoper::CandidateFinder(NumericVector PeaksToCheck)
{
  double* limits; //Left and rigth mass limits for each candidate
  int candidateCnt = 0;    //counter of candidates for each peak
  int candidateIndx = 0;  //Index of the candidate peak
  
  //Matrix containing the candidate masses for each selected peak
  NumericMatrix CanMatrix(imgRuninfo->massPeaks, 21); //Allowed a maximum of 20 candidates
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    for(int j = 0; j < 21; j++)
    {
      CanMatrix(i,j) = 0; 
    }
  }
  
  //Filling the first column of the candidates matrix with the peak masses index
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    CanMatrix(i,0) = pMatrixPeakOrder[i];
  }
  
  //Searching the candidates peak indexes of each mass
  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    if(PeaksToCheck[i] > 0)
    {
      limits = getCandidateLimits(i);
      candidateCnt = 1; 
      for(int j = 1; j < imgRuninfo->massPeaks; j++)
      {
        candidateIndx = pMatrixPeakOrder[i] + j;
        if(candidateIndx >= imgRuninfo->massPeaks)    //Outside matrix
        {
          break;
        }
        
        if( ((pMatrixAxis[candidateIndx]) > (limits[0])) &&   //Between limits
            ((pMatrixAxis[candidateIndx]) < (limits[1])) ) 
        {
          CanMatrix(i, candidateCnt) = candidateIndx;
          candidateCnt++;
        }
        
        if((pMatrixAxis[candidateIndx]) > (limits[1]))  //Outside limits & no more possible candidates
        {
          break;
        }
      }
      pNumCandidates[i] = candidateCnt-1; //We start to cont from 1 so if 0 candidates 1-1. just to not overwrite the can matrix 0 column
    }
  }
  
  //Finding the maximum nubmer of candidates a mass has
  maxCan = 0;

  for(int i = 0; i < imgRuninfo->massPeaks; i++)
  {
    if(maxCan < pNumCandidates[i])
    {
      maxCan = pNumCandidates[i];
    }
  }

  return CanMatrix;
}

double* Deisotoper::getCandidateLimits(int massIndex)
{
  static double limits[2] = {0,0};
  double CentralIsoMass = pMatrixAxis[pMatrixPeakOrder[massIndex]] + (CrntNumIso*1.0033548378);
  double tmpDummy = 0; //If scans it helps to find the closest scan, if ppm it stores the ppm error from the central mass
  
  if(imgRuninfo->ToleranceInScans)
  {
    int ImageAxisPointer = 0;
    
    tmpDummy = pImageAxis[imgRuninfo->massChannels - 1] - pImageAxis[0]; //Maximum distance between scans
    
    for(int k = pImagePeakOrder[massIndex]; k < imgRuninfo->massChannels; k++) //Founds the scan closer to CentralIsoMass and stores it in ImageAxisPointer
    {
      if(std::abs(CentralIsoMass - pImageAxis[k]) > tmpDummy)
      {
        ImageAxisPointer = k - 1;
        break;
      }
      
      if(std::abs(CentralIsoMass - pImageAxis[k]) < tmpDummy)
      {
        tmpDummy = std::abs(CentralIsoMass - pImageAxis[k]);
      }
    }
    
    limits[0] = ((ImageAxisPointer - imgRuninfo->tolerance) < 0) ? pImageAxis[0] : pImageAxis[ImageAxisPointer - imgRuninfo->tolerance];
    limits[1] = ((ImageAxisPointer + imgRuninfo->tolerance) > imgRuninfo->massChannels) ? pImageAxis[imgRuninfo->massChannels] : pImageAxis[ImageAxisPointer + imgRuninfo->tolerance];
  } 
    else
    {
      tmpDummy = imgRuninfo->tolerance*CentralIsoMass;
      tmpDummy = tmpDummy/1000000;
      limits[0] = CentralIsoMass - tmpDummy;  //Left limit
      limits[1] = CentralIsoMass + tmpDummy;  //Right limit
    }
  
  return(limits);
}

//Gets the tolerance from the mass error curve calculated during the constructor
double Deisotoper::getToleranceFromCurve(int imgIndex)
{
  return(massErrorCurve[imgIndex]);
}

//Computes all the scores over a candidates row 
double* Deisotoper::ScoreCalculator(int* CandidateRow, int NumCan, double* result, double lastModslope, int peakNumber)
{
  double x[imgRuninfo->numPixels];  //array containing the image of the principal mass 
  double y[imgRuninfo->numPixels];  //array containing the image of the candidate mass
  double A = 0, B = 0, y_mean = 0, x_mean = 0, SStot, SSres, Modslope, intercept;
  double ScoreMrph, ScoreInt,ScoreMass, ModCA, CA = 0,ppm = 0 ,maxppm;
  
  for(int i = 0; i < (7*NumCan); i++) //Cleanning 
  {
    result[i] = 0;
  }
  
  for(int i = 0; i < imgRuninfo->numPixels ; i++)
  {
    x[i] = pMatrix[i][CandidateRow[0]];  //Monoisotopic peak
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
      y[j] = pMatrix[j][CandidateRow[i+1]];  //M+N peak
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
      Modslope = (A/B < 0) ?  0 : A/B;  //Isotope ratio from the data
      
      intercept = y_mean - Modslope*x_mean;
      
    //***********************************//Morphology score (adj. R2)//*******************************************************//
    
      //Total sum of squares: the sum of the squares of the difference of the dependent variable and its mean
      for(int j = 0; j < imgRuninfo->numPixels ; j++)
      {
        SStot += (y[j] - y_mean) * (y[j] - y_mean);     
      }
      
      //Sum of squared residuals: the sum of the squares of the residuals
      for(int j = 0; j < imgRuninfo->numPixels ; j++)
      {
        SSres += (y[j] - (intercept + (Modslope*x[j]))) * (y[j] - (intercept + (Modslope*x[j])));     //TODO probar treure el model i ficar directament el valor de x
      }
      
      //Definition of R2 
      ScoreMrph = 1 - (SSres/SStot); 
      ScoreMrph = sqrt(ScoreMrph);
      
    //***********************************//Intensity Score//******************************************************************//
    
      //Theoretical number of Carbons extracted from the Human metabolome database model  
      ModCA  = 0.076*pMatrixAxis[CandidateRow[0]] - 12;
      
      //Number of carbons from the experimental intensity rate
      
      CA = (CrntNumIso == 2) ? ((sqrt(1+(8*Modslope*(0.9893/0.0107)*(0.9893/0.0107)))+1)*0.5) : (CrntNumIso*((Modslope/lastModslope)*(0.9893/0.0107) + 1) - 1);

      ScoreInt  = abs(ModCA - CA)/(sqrt(1/(-log(0.7)))*20);  //Adjusting 20 atoms of difference to score 0.7
      ScoreInt *= -ScoreInt;                    
      ScoreInt  = exp(ScoreInt); 
    
    //***********************************//Mass ppm Score//******************************************************************//
    
    ppm = abs((((pMatrixAxis[CandidateRow[i+1]])-(pMatrixAxis[CandidateRow[0]]+CrntNumIso*1.0033548378))*1000000)/(pMatrixAxis[CandidateRow[0]]+CrntNumIso*1.0033548378)); //mass error in ppm
    
    maxppm = (!imgRuninfo->ToleranceInScans) ? imgRuninfo->tolerance : getToleranceFromCurve(pImagePeakOrder[peakNumber]); //maximum valid tolerance
  
    ScoreMass  =  ppm/(sqrt((1/log(2)))*maxppm); //Adjusting the maximum tolerance to score 0.5
    ScoreMass *= -ScoreMass;
    ScoreMass  = exp(ScoreMass);
    
    //***********************************//Results//*******************************************************//  
    
      //Returning the scores
      result[(7*i) + 0] = ppm;                                    //ppm error
      result[(7*i) + 1] = ((ScoreMrph + ScoreInt + ScoreMass)/3); //final score 
      result[(7*i) + 2] = ScoreMrph;                              //morphology score
      result[(7*i) + 3] = ScoreInt;                               //intensity score
      result[(7*i) + 4] = ScoreMass;                              //mass error score
      result[(7*i) + 5] = Modslope;                                //model slope  
      result[(7*i) + 6] = CA;                                     //number of C atmos
  }

  return result;
}

//Directs the ScoreCalculator function for each ion to its candidates
List Deisotoper::MatrixAnnotator(NumericMatrix CanMatrix, NumericVector PeaksToCheck)
{
  List AnnMatrix(imgRuninfo->massPeaks);   // List containing the results of the annotation process
  
  int* CandidateRow;  //pointer to a vector containing the index of the masses to be tested
  CandidateRow = new int[maxCan];
  
  double* scores;   //pointer to a vector containing the results of the ScoreCalculator function
  scores = new double[7*maxCan];
  for(int i = 0; i < 7*maxCan; i++)
  {
    scores[i] = 0;
  }
  
  NumericMatrix TMP(10,maxCan); 
  
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
  List Result(imgRuninfo->numIso + 1);  //List containing the output results
  List PeakResults; //Results for each isotopoic stage
  NumericMatrix PeaksToCheck(imgRuninfo->massPeaks, imgRuninfo->numIso); //Matrix that indicates if a ion must undergo the test
                                                                         //A 0 means not to check
                                                                         //This matrix is overwritten with the slope of the previous stage if the score is higher than score threshold
  //***********************************//Fill the PeaksToCheck matrix//***********************************************// 
  for(int j = 0; j < imgRuninfo->numIso; j++)  
  {
    if(j == 0)
      {
        for(int i = 0; i < imgRuninfo->massPeaks; i++)
        {
          PeaksToCheck(i,j) = 1;  //For the first stage (M1), all peaks are gonna be checked
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
  double tHighestScore = 0; //Highest score at that time
  for(int i = 0; i < imgRuninfo->numIso; i++)  
  {
    CrntNumIso = i + 1;  //Current searched isotope number
    NumericVector PksTChk = PeaksToCheck(_,i); //Current PeaksToCheck, contains information of the previous stage.

    //Matrix filled with the candidates for each monoisotopic peak
    NumericMatrix CandidatesMatrix = CandidateFinder(PksTChk);  //Generates a matrix that contains all the candidats to be isotopes for each peak mass
    PeakResults = MatrixAnnotator(CandidatesMatrix, PksTChk); //Computes the scores for each candidate in the candidate matrix that are 1st isotopes or have a positive score in the previous stage.
    Result[i + 1] = PeakResults; //Results[0] containts the sorted peak matrix mass aixs
    
    if((CrntNumIso < imgRuninfo->numIso)) //If there are more stages, fill the PeaksToCheck
    {
      for(int j = 0; j < imgRuninfo->massPeaks; j++) //Filling the PeksToCheck vector with the new results
      {
        if((PksTChk[j] > 0) & (pNumCandidates[j] > 0)) //Check if the tests has been done and if there are any candidates
        {
          NumericMatrix TestResults = PeakResults[j];
          tHighestScore = 0;
          for(int k = 0; k < pNumCandidates[j]; k++)
          {
            if((TestResults(4,k) > imgRuninfo->scoreThreshold) & (TestResults(4,k) > tHighestScore)) //Checks if the candidate has passed the test & if its the candidate with higher score
            {
              PeaksToCheck(j,CrntNumIso) = TestResults(8,k); //Slope of the previous step 
              tHighestScore = TestResults(8,k);
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
Rcpp::List C_isotopeAnnotator(int massPeaks, int massChannels, int numPixels, int numIso,
                              NumericMatrix PeakMtx, NumericVector massVec, NumericVector massChanVec,
                              int tolerance, double scoreThreshold, bool ToleranceInScans)
{
  //Fill the class data structure with the information of the experiment
  Deisotoper::IsoDef myIsoDef;  
  myIsoDef.massPeaks = massPeaks;
  myIsoDef.massChannels = massChannels;
  myIsoDef.numPixels = numPixels;
  myIsoDef.numIso = numIso;
  myIsoDef.tolerance = tolerance;
  myIsoDef.scoreThreshold = scoreThreshold;
  myIsoDef.ToleranceInScans = ToleranceInScans;
  
  List result;
  Deisotoper myDeisotoper(&myIsoDef, PeakMtx, massVec, massChanVec); //Class constructor
  result = myDeisotoper.Run();  //Algorithm run
  
  return result;
}

