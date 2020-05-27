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

#include "adductannotation.h"
using namespace Rcpp;

//Constructor
AdductPairer::AdductPairer(AdductDef *imgRunInfo, NumericVector R_monoisitopeMassVector,
                           NumericVector R_adductMassVector, List R_isotopes,
                           NumericVector R_isotopeListOrder, NumericVector R_massAxis,
                           NumericMatrix R_peakMatrix, NumericVector R_labelAxis,
                           NumericVector R_monoisotopicIndexVector)
{
	RunDef = imgRunInfo;
  positiveTest = 0;
  isotopes = R_isotopes;
  
  
	monoisitopeMassVector = new double[RunDef->numMonoiso];
	for(int i = 0; i < RunDef->numMonoiso; i++)
	{
		monoisitopeMassVector[i] = R_monoisitopeMassVector[i];
	}
	
	monoisitopeIndexVector = new int[RunDef->numMonoiso];
	for(int i = 0; i < RunDef->numMonoiso; i++)
	{
	  monoisitopeIndexVector[i] = R_monoisotopicIndexVector[i];
	}

	
	massAxisToMonoisotopicAxis = new int[RunDef->numMass];
	for(int i = 0; i < RunDef->numMass; i++)
	{
	  massAxisToMonoisotopicAxis[i] = -1;
	}
	int cnt = 0;
	for(int i = 0; i < RunDef->numMonoiso; i++)
	{
	  massAxisToMonoisotopicAxis[monoisitopeIndexVector[i]] = cnt;
	  cnt++;
	}
	  
	massAxis = new double[RunDef->numMass];
	for(int i = 0; i < RunDef->numMass; i++)
	{
	  massAxis[i] = R_massAxis[i];
	}
	
	
	adductMassVector = new double[RunDef->numAdducts];
	for(int i = 0; i < RunDef->numAdducts; i++)
	{
	  adductMassVector[i] = R_adductMassVector[i];
	}
	
	
	isotopeListOrder = new int[RunDef->numMonoiso];
	for(int i = 0; i < RunDef->numMonoiso; i++)
	{
	  isotopeListOrder[i] = R_isotopeListOrder[i];
	}
	
	
	labelAxis = new int[RunDef->numMass];
	for(int i = 0; i < RunDef->numMass; i++)
	{
	  labelAxis[i] = R_labelAxis[i];
	}
  
  peakMatrix = new double*[RunDef->numPixels];
  for(int i = 0; i < RunDef->numPixels; i++)
  {
    peakMatrix[i] = new double[RunDef->numMass];
    for(int j = 0; j < RunDef->numMass; j++)
    {
      peakMatrix[i][j] = R_peakMatrix(i,j);
    }
  }
  
}


//Destructor
AdductPairer::~AdductPairer()
{
	delete[] monoisitopeMassVector;
  
  delete[] monoisitopeIndexVector;
  
  delete[] massAxisToMonoisotopicAxis;
  
  delete[] labelAxis;
  
  for(int i = 0; i < RunDef->numMass; i++)
  {
    delete[] ppmMatrix[i];
  }
  delete[] ppmMatrix;
  
  
  for(int i = 0; i < RunDef->numPixels; i++)
  {
    delete[] peakMatrix[i];
  }
  delete[] peakMatrix;
  
  for(int i = 0; i < RunDef->numMass; i++)
  {
    delete[] neutralMassMatrix[i];
  }
  delete[] neutralMassMatrix;
  
  for(int i = 0; i < RunDef->numMass; i++)
  {
    delete[] labelMatrix[i];
  }
  delete[] labelMatrix;
  
	delete[] adductMassVector;
	
	
	delete[] isotopeListOrder;
	
	
	delete[] massAxis;
	

	delete[] adductPairFirstMassVector;
	
	
	delete[] adductPairSecondMassVector;
}


//Creats the two vectors were the masses of the elements are stored and creates the ppmmatrix, massMatrix and labelmatrix
void AdductPairer::AdductPairCreator()
{
	adductDiffCombinations = 0;	  //Calculating the number of pair of adducts combinations
	for(int i = 0; i < RunDef->numAdducts; i++)
	{
		adductDiffCombinations += i;
	}

  adductPairFirstMassVector = new double[adductDiffCombinations];
  adductPairSecondMassVector = new double[adductDiffCombinations];
  
	int p1 = 0;   //dummy counter
	int p2 = 1;   //dummy counter
	int lim = RunDef->numAdducts-1;   //dummy flag

	for(int i = 0; i < adductDiffCombinations; i++)   //Fill the vector
	{
	  adductPairFirstMassVector[i] = adductMassVector[p1];
	  adductPairSecondMassVector[i] = adductMassVector[p2];
		p2++;
		if(p2 >lim)
		{
			p1++;
			p2 = p1+1;
		}
	}
	
	neutralMassMatrix = new double*[RunDef->numMass];
	for(int i = 0; i < RunDef->numMass; i++)
	{
	  neutralMassMatrix[i] = new double[adductDiffCombinations*RunDef->numMonoiso];
	  for(int j = 0; j < (adductDiffCombinations*RunDef->numMonoiso); j++)
	  {
	    neutralMassMatrix[i][j] = 0;
	  }
	}
	
	ppmMatrix = new double*[RunDef->numMass];
	for(int i = 0; i < RunDef->numMass; i++)
	{
	  ppmMatrix[i] = new double[adductDiffCombinations*RunDef->numMonoiso];
	  for(int j = 0; j < (adductDiffCombinations*RunDef->numMonoiso); j++)
	  {
	    ppmMatrix[i][j] = 9999;
	  }
	}
	
	labelMatrix = new int*[RunDef->numMass];
	for(int i = 0; i < RunDef->numMass; i++)
	{
	  labelMatrix[i] = new int[adductDiffCombinations*RunDef->numMonoiso];
	  for(int j = 0; j < (adductDiffCombinations*RunDef->numMonoiso); j++)
	  {
	    labelMatrix[i][j] = -1;
	  }
	}
}



//Substracts the monoisotopeMassVector with the shifted version of himself
void AdductPairer::CalculateAdductPairs()
{
  double ppmMassError = 9999; //ppm mass error for each pair
  double Mmass1;
  double Mmass2;
  double meanMmass;
  double mass_diff;

    for(int i = 0; i < adductDiffCombinations; i++) 
    {
      for(int j = 0; j < RunDef->numMonoiso; j++) 
      {
        for(int k = 0; k < RunDef->numMass; k++)
        {
          if(labelAxis[k] != 2)   //check if the peak is not an isotopic peak
          {
            if(monoisitopeMassVector[j] > massAxis[k])
            {
              Mmass1 = monoisitopeMassVector[j] - adductPairFirstMassVector[i];
              Mmass2 = massAxis[k] - adductPairSecondMassVector[i];
              meanMmass = (Mmass1 + Mmass2)/2;
              mass_diff = Mmass1 - Mmass2; 
              ppmMassError  = 1e6*mass_diff/meanMmass;
              if(fabs(ppmMassError) < RunDef->tolerance)
              {
                ppmMatrix[k][(i*RunDef->numMonoiso)+j] = fabs(ppmMassError);
                neutralMassMatrix[k][(i*RunDef->numMonoiso)+j] = meanMmass;
                labelMatrix[k][(i*RunDef->numMonoiso)+j] = labelAxis[k];
                positiveTest++;
              }
            }
            //else if((monoisitopeMassVector[j] < massAxis[k]))
            else if((monoisitopeMassVector[j] < massAxis[k]) & (labelAxis[k] == 0))
            {
              Mmass1 = monoisitopeMassVector[j] - adductPairSecondMassVector[i];
              Mmass2 = massAxis[k] - adductPairFirstMassVector[i];
              meanMmass = (Mmass1 + Mmass2)/2;
              mass_diff = Mmass1 - Mmass2; 
              ppmMassError  = 1e6*mass_diff/meanMmass;
              if(fabs(ppmMassError) < RunDef->tolerance)
              {
                ppmMatrix[k][(i*RunDef->numMonoiso)+j] = fabs(ppmMassError);
                neutralMassMatrix[k][(i*RunDef->numMonoiso)+j] = meanMmass;
                labelMatrix[k][(i*RunDef->numMonoiso)+j] = labelAxis[k];
                positiveTest++;
              }
            }
          }
        }
      }
    }  
}


//Reads the information form the isotopes tests and merge it with the adduct candidates information.
// void AdductPairer::ValidateCandidates()
// {
//   double slope1 = 0;
//   double slope2 = 0;
//   List IsoData1;
//   List IsoData2;
//   int mdVlength;
//   
//   for(int i = 0; i < maxShift; i++)
//   {
//     mdVlength = RunDef->numMonoiso - 1 - i;
//     for(int j = 0; j < mdVlength; j++)
//     {
//       if(adductPairsNameMatrix[i][j] != -1)
//       {
//         IsoData1 = isotopes[isotopeListOrder[j+i+1]];
//         IsoData2 = isotopes[isotopeListOrder[j]];
//         slope1 = IsoData1[8];
//         slope2 = IsoData2[8];
//         
//         if(abs(slope1 - slope2) > 0.2)    //positive test condition (may change)
//         {
//           positiveTest -= 1;
//           adductPairsNameMatrix[i][j] = -1;
//         } 
//       }
//     }
//   }
// }


//Output construction function
List AdductPairer::GenerateResultsList()
{
  int cntResultsA = 0;
  int cntResultsB = 0;
  List Results(2);
	List ResultsA(positiveTest+1);
	List ResultsB(positiveTest+1);
	//NumericVector names(positiveTest);
	NumericVector peakA(RunDef->numPixels), peakB(RunDef->numPixels); 
	double tmp = 0;
	for(int i = 0; i < adductDiffCombinations; i++) 
	{
	  for(int j = 0; j < RunDef->numMonoiso; j++) 
	  {
	    for(int k = 0; k < RunDef->numMass; k++)
	    {
  		  if(labelMatrix[k][(i*RunDef->numMonoiso)+j] == 1)
  		  {
  		    List IsoData1 = isotopes(isotopeListOrder[j]);
  		    List IsoData2 = isotopes(isotopeListOrder[massAxisToMonoisotopicAxis[k]]);
  		    NumericVector r(10);

  		    r(0) = i;                               //Pair of adducts ID
  		    r(1) = neutralMassMatrix[k][(i*RunDef->numMonoiso)+j]; //Neutral mass
  		    r(2) = monoisitopeMassVector[j];        //M+Add1 mass
  		    r(3) = IsoData1(2);                     //M+Add1 index
  		    r(4) = monoisitopeMassVector[massAxisToMonoisotopicAxis[k]];        //M+Add1 mass
  		    r(5) = IsoData2(2);                     //M+Add1 index
  		    
  		    tmp = IsoData1(9);                        //Old: Average C atoms. New: Average slopes 
  		    tmp *= (0.0107/0.9893);
  		    r(6)  = IsoData2(9);  
  		    r(6) *= (0.0107/0.9893);
  		    r(6) += tmp;   
  		    r(6) *= 0.5;  
  		    
  	      tmp = IsoData1(9);                        //C atmos error * probabilities quotien = M+1 slope
  	      tmp *= (0.0107/0.9893);
  	      tmp = (tmp-r(6))*(tmp-r(6));
  		    r(7) = IsoData2(9); 
  		    r(7) *= (0.0107/0.9893);
  		    r(7) = (r(7)-r(6))*(r(7)-r(6));
  		    r(7) = sqrt((r(7) + tmp))/sqrt(2);
  		    
  		    for(int l = 0; l < RunDef->numPixels; l++)
  		    {
  		      peakA[l] = peakMatrix[l][(int)r(3)-1];
  		      peakB[l] = peakMatrix[l][(int)r(5)-1];
  		    }
  
  		    r(8) = correlation(peakA, peakB, RunDef->numPixels);   //Correlation degree
  		    r(9) = ppmMatrix[k][(i*RunDef->numMonoiso)+j];         //Mass diference error in ppm   
  		    
  		    ResultsA(cntResultsA) = r;
  		    cntResultsA++;
  		    
  		    // names[cntList]  = massMatrix[i][j];
		    }
  		  
  		  if(labelMatrix[k][(i*RunDef->numMonoiso)+j] == 0)
  		  {
  		    List IsoData1 = isotopes(isotopeListOrder[j]);
  		    NumericVector r(8);
  		    r(0) = i;                                                 //Adduct pair code
  		    r(1) = neutralMassMatrix[k][(i*RunDef->numMonoiso)+j];    //Neutral mass
  		    r(2) = monoisitopeMassVector[j];                          //Monoisotopic mass
  		    r(3) = IsoData1(2);                                       //Monoisotopic mass index THIS SHOULD BE WRONG
  		    r(4) = massAxis[k];                                       //Unknown mass
  		    r(5) = k+1;                                               //Unknown mass index
  		    
  		    for(int l = 0; l < RunDef->numPixels; l++)
  		    {
  		      peakA(l) = peakMatrix[l][(int)r(3)-1];
  		      peakB(l) = peakMatrix[l][(int)r(5)-1];
  		    }
  	
  		    r(6) = correlation(peakA,peakB,RunDef->numPixels);   //Correlation degree
  		    r(7) = ppmMatrix[k][(i*RunDef->numMonoiso)+j];                        //Mass diference error in ppm   
  		    
  		    ResultsB(cntResultsB) = r;
  		    cntResultsB++;
  		  }
  		  
  		  if(cntResultsA+cntResultsB == positiveTest)
  		  {
  		    Results(0) = ResultsA;
  		    Results(1) = ResultsB;
  		    return(Results);
  		  }
		  }
		}
	}

	return(Results);
}


//Program executioning function
List AdductPairer::Run()
{
  AdductPairCreator();                //Fills adductMassDifferencesVector

  CalculateAdductPairs();            //Do the shift and substract the masses

  //ValidateCandidates();           //Reads the information form the isotopes tests and merge it with the adduct candidates information.
	
	Rcout << "Number of adduct pairs found = "<< positiveTest << "\n";
	
	return(GenerateResultsList());
}                                                                   


// [[Rcpp::export]]
Rcpp::List C_adductAnnotation(int numMonoiso, int numAdducts, 
                              int tolerance, int numMass, NumericVector R_monoisitopeMassVector, 
                              NumericVector R_adductMassVector, List R_isotopes,
                              NumericVector R_isotopeListOrder, NumericVector R_massAxis, 
                              NumericMatrix R_peakMatrix, int numPixels, NumericVector R_labelAxis,
                              NumericVector R_monoisotopicIndexVector)
{
  //Fill the class data structure with the information of the experiment
  AdductPairer::AdductDef	myAdductDef;  
  myAdductDef.numMass = numMass;
  myAdductDef.numAdducts = numAdducts;
  myAdductDef.numMonoiso = numMonoiso;
  myAdductDef.tolerance = tolerance;
  myAdductDef.numPixels = numPixels;

  List result;
  AdductPairer myAdductPairer(&myAdductDef, R_monoisitopeMassVector, R_adductMassVector, R_isotopes, R_isotopeListOrder, R_massAxis, R_peakMatrix, R_labelAxis, R_monoisotopicIndexVector); //Class constructor
  
  result = myAdductPairer.Run();  //Program executioning function
  
  return result;
}


