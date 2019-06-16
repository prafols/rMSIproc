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
                           NumericMatrix R_peakMatrix)
{
	RunDef = imgRunInfo;
  currentShift = 0;
  positiveTest = 0;
  mdVlength = RunDef->numMonoiso - 1;
  maxShift = RunDef->numMonoiso - 1;
  isotopes = R_isotopes;
  
  
	monoisitopeMassVector = new double[RunDef->numMonoiso];
	for(int i = 0; i < RunDef->numMonoiso; i++)
	{
		monoisitopeMassVector[i] = R_monoisitopeMassVector[i];
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
	
	
	ppmMatrix = new double*[maxShift];
  for(int i = 0; i < maxShift; i++)
  {
    ppmMatrix[i] = new double[mdVlength];
    for(int j = 0; j < mdVlength; j++)
    {
      ppmMatrix[i][j] = 9999;
    }
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
  
  massMatrix = new double*[maxShift];
  for(int i = 0; i < maxShift; i++)
  {
    massMatrix[i] = new double[mdVlength];
    for(int j = 0; j < mdVlength; j++)
    {
      massMatrix[i][j] = 0;
    }
  }

  adductPairsNameMatrix = new int*[maxShift];
  for(int i = 0; i < maxShift; i++)
  {
    adductPairsNameMatrix[i] = new int[mdVlength];
    for(int j = 0; j < mdVlength; j++)
    {
      adductPairsNameMatrix[i][j] = -1;	//-1 means that the combination is not related with adducts of any type
    }
  }
}


//Destructor
AdductPairer::~AdductPairer()
{
	delete[] monoisitopeMassVector;
  
  
  for(int i = 0; i < maxShift; i++)
  {
    delete[] ppmMatrix[i];
  }
  delete[] ppmMatrix;
  
  
  for(int i = 0; i < maxShift; i++)
  {
    delete[] adductPairsNameMatrix[i];
  }
  delete[] adductPairsNameMatrix;
  
  for(int i = 0; i < RunDef->numPixels; i++)
  {
    delete[] peakMatrix[i];
  }
  delete[] peakMatrix;
  
  for(int i = 0; i < maxShift; i++)
  {
    delete[] massMatrix[i];
  }
  delete[] massMatrix;
  
	delete[] adductMassVector;
	
	
	delete[] isotopeListOrder;
	
	
	delete[] massAxis;
	
	
	delete[] adductPairFirstMassVector;
}


//Fills adductMassDifferencesVector
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
}


//Substracts the monoisotopeMassVector with the shifted version of himself
void AdductPairer::ShiftandCalculateAdductPairs()
{
  double ppmMassError = 9999; //ppm mass error for each pair
  double Mmass1;
  double Mmass2;
  double meanMmass;
  double mass_diff;
  currentShift = currentShift + 1;  //Add more shift
  mdVlength = RunDef->numMonoiso - currentShift;  //Shrink the massdifferences matrix row
  
    for(int i = 0; i < adductDiffCombinations; i++) 
    {
      for(int j = 0; j < mdVlength; j++)  //Fill the row with the differences
      {
        Mmass1 = monoisitopeMassVector[currentShift+j] - adductPairFirstMassVector[i];
        Mmass2 = monoisitopeMassVector[j] - adductPairSecondMassVector[i];
        meanMmass = (Mmass1 + Mmass2)/2;
        mass_diff = Mmass1 - Mmass2; 
        ppmMassError  = 1e6*mass_diff/meanMmass;
        if(fabs(ppmMassError) < RunDef->tolerance)
        {
          adductPairsNameMatrix[currentShift-1][j] = i;
          ppmMatrix[currentShift-1][j] = ppmMassError;
          massMatrix[currentShift-1][j] = meanMmass;
          positiveTest++;
        }
      }
    }  
}


//Reads the information form the isotopes tests and merge it with the adduct candidates information.
void AdductPairer::ValidateCandidates()
{
  double slope1 = 0;
  double slope2 = 0;
  List IsoData1;
  List IsoData2;
  int mdVlength;
  
  for(int i = 0; i < maxShift; i++)
  {
    mdVlength = RunDef->numMonoiso - 1 - i;
    for(int j = 0; j < mdVlength; j++)
    {
      if(adductPairsNameMatrix[i][j] != -1)
      {
        IsoData1 = isotopes[isotopeListOrder[j+i+1]];
        IsoData2 = isotopes[isotopeListOrder[j]];
        slope1 = IsoData1[8];
        slope2 = IsoData2[8];
        
        if(abs(slope1 - slope2) > 0.2)    //positive test condition (may change)
        {
          positiveTest -= 1;
          adductPairsNameMatrix[i][j] = -1;
        } 
      }
    }
  }
}


//Output construction function
List AdductPairer::GenerateResultsList()
{
  int cntList = 0;
	List Results(positiveTest+1);
	NumericVector names(positiveTest);
	NumericVector peakA(RunDef->numPixels), peakB(RunDef->numPixels); 
	double tmp = 0;
	
	for(int i = 0; i < maxShift; i++)
	{
	  mdVlength = RunDef->numMonoiso - i - 1;
		for(int j = 0; j < mdVlength; j++)
		{
		  if (adductPairsNameMatrix[i][j] != -1)
		  {
		    List IsoData1 = isotopes[isotopeListOrder[j+(i+1)]];
		    List IsoData2 = isotopes[isotopeListOrder[j]];
		    NumericMatrix r(9,1);
		    
		    r(0,0) = adductPairsNameMatrix[i][j];     //Pair of adducts ID
		    r(1,0) = monoisitopeMassVector[j+i+1];    //M+Add1 mass
		    r(2,0) = IsoData1[2];                     //M+Add1 index
		    r(3,0) = monoisitopeMassVector[j];        //M+Add1 mass
		    r(4,0) = IsoData2[2];                     //M+Add1 index
		    
		    tmp = IsoData1[9];                        //Average C atoms
		    r(5,0)  = IsoData2[9];                    
		    r(5,0) += tmp;   
		    r(5,0) *= 0.5;  
	
	      tmp = IsoData1[9];                        //C atmos error
	      tmp = abs(tmp-r(5,0));
		    r(6,0) = IsoData2[9];                       
		    r(6,0) = abs(r(6,0)-r(5,0));
		    r(6,0) = (r(6,0) + tmp)*0.5;
		    
		    for(int k = 0; k < RunDef->numPixels; k++)
		    {
		      peakA[k] = peakMatrix[k][(int)r(2,0)-1];
		      peakB[k] = peakMatrix[k][(int)r(4,0)-1];
		    }

		    r(7,0) = correlation(peakA,peakB,RunDef->numPixels);   //Correlation degree
		    r(8,0) = ppmMatrix[i][j];                              //Mass diference error in ppm   
		    
		    
		    Results[cntList] = r;

		    names[cntList]  = massMatrix[i][j];
		    cntList = cntList + 1;
		  }
		}
	}
	Results[cntList] = names;
	return(Results);
}


//Program executioning function
List AdductPairer::Run()
{
  AdductPairCreator(); //Fills adductMassDifferencesVector

	while(currentShift<maxShift)      //While more shifts allowed
	{
	  ShiftandCalculateAdductPairs();            //Do the shift and substract the masses
	}
	
  //ValidateCandidates();           //Reads the information form the isotopes tests and merge it with the adduct candidates information.
	
	Rcout << "Number of adduct pairs found = "<< positiveTest << "\n";
	
	return(GenerateResultsList());
}                                                                   


// [[Rcpp::export]]
Rcpp::List C_adductAnnotation(int numMonoiso, int numAdducts, 
                              int tolerance, int numMass, NumericVector R_monoisitopeMassVector, 
                              NumericVector R_adductMassVector, List R_isotopes,
                              NumericVector R_isotopeListOrder, NumericVector R_massAxis, 
                              NumericMatrix R_peakMatrix, int numPixels)
{
  //Fill the class data structure with the information of the experiment
  AdductPairer::AdductDef	myAdductDef;  
  myAdductDef.numMass = numMass;
  myAdductDef.numAdducts = numAdducts;
  myAdductDef.numMonoiso = numMonoiso;
  myAdductDef.tolerance = tolerance;
  myAdductDef.numPixels = numPixels;
  
  List result;
  AdductPairer myAdductPairer(&myAdductDef, R_monoisitopeMassVector, R_adductMassVector, R_isotopes, R_isotopeListOrder, R_massAxis, R_peakMatrix); //Class constructor
  
  result = myAdductPairer.Run();  //Program executioning function
  
  return result;
}


