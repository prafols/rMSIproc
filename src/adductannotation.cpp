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
AdductPairer::AdductPairer(AdductDef *imgRunInfo, NumericVector R_monoisitopeMassVector, NumericVector R_adductMassVector, List R_isotopes, NumericVector R_isotopeListOrder, NumericVector R_massAxis)
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
  
  
	delete[] adductMassVector;
	
	
	delete[] isotopeListOrder;
	
	
	delete[] massAxis;
	
	
	delete[] adductPairFirstMassVector;
}


//Fills adductMassDifferencesVector
void AdductPairer::CalculateAdductMassDifferences()
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
  double ppmMassError = 0; //ppm mass error for each pair
  double Mmass1;
  double Mmass2;
  
  currentShift = currentShift + 1;  //Add more shift
  mdVlength = RunDef->numMonoiso - currentShift;  //Shrink the massdifferences matrix row
  
  for(int i = 0; i < adductDiffCombinations; i++) 
  {
    for(int j = 0; j < mdVlength; j++)  //Fill the row with the differences
    {
      Mmass1 = monoisitopeMassVector[currentShift+j]-adductPairFirstMassVector[i];
      Mmass2 = monoisitopeMassVector[j]-adductPairSecondMassVector[i];
      ppmMassError  = 1000000*abs(Mmass1 - Mmass2)/(Mmass2);
      adductPairsNameMatrix[currentShift-1][j] = ((ppmMassError < RunDef->tolerance) & (adductPairsNameMatrix[currentShift-1][j] == -1)) ? i : -1;  
      ppmMatrix[currentShift-1][j] = (adductPairsNameMatrix[currentShift-1][j] != -1) ? ppmMassError : ppmMatrix[currentShift-1][j];  
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
  
  for(int i = 0; i < maxShift; i++)
  {
    mdVlength = RunDef->numMonoiso - (i+1);
    for(int j = 0; j < mdVlength; j++)
    {
      if(adductPairsNameMatrix[i][j] != -1)
      {
        IsoData1 = isotopes[isotopeListOrder[j+(i+1)]];
        IsoData2 = isotopes[isotopeListOrder[j]];
        slope1 = IsoData1[8];
        slope2 = IsoData2[8];
        
        if(abs(slope1 - slope2) < 0.3)    //positive test condition (may change)
        {
          positiveTest += 1;
        } 
        else
        {
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
	  
	for(int i = 0; i < maxShift; i++)
	{
		for(int j = 0; j < (RunDef->numMonoiso-1); j++)
		{
		  if (adductPairsNameMatrix[i][j] != -1)
		  {
		    List IsoData1 = isotopes[isotopeListOrder[j+(i+1)]];
		    List IsoData2 = isotopes[isotopeListOrder[j]];
		    NumericMatrix r(9,1);
		    
		    r(0,0)  =  adductPairsNameMatrix[i][j];         //Pair of adducts ID
		    
		    r(1,0)  = massAxis[as<int>(IsoData1[2])-1];     //M+Add1 mass
		    r(2,0)  = IsoData1[2];                          //M+Add1 index
		    r(3,0)  = IsoData1[8];                          //M+Add1 slope
		    
		    r(4,0)  = massAxis[as<int>(IsoData2[2])-1];     //M+Add1 mass
		    r(5,0)  = IsoData2[2];                          //M+Add1 index
		    r(6,0)  = IsoData2[8];                          //M+Add2 slope
		    
		    r(7,0)  = j;                             //Correlation degree
		    r(8,0)  = ppmMatrix[i][j];   //Mass diference error in ppm 
		    
		    Results[cntList] = r;

		    names[cntList]  = (massAxis[as<int>(IsoData1[2])-1] - adductPairFirstMassVector[adductPairsNameMatrix[i][j]]) ;    //Absolute mass of the compound
		    names[cntList] += (massAxis[as<int>(IsoData2[2])-1] - adductPairSecondMassVector[adductPairsNameMatrix[i][j]]);
		    names[cntList] /= 2;
		    
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
	CalculateAdductMassDifferences(); //Fills adductMassDifferencesVector

	while(currentShift<maxShift)      //While more shifts allowed
	{
	  ShiftandCalculateAdductPairs();            //Do the shift and substract the masses
	}


	ValidateCandidates();           //Reads the information form the isotopes tests and merge it with the adduct candidates information.
	
	return(GenerateResultsList());
}                                                                   


// [[Rcpp::export]]
Rcpp::List C_adductAnnotation(int numMonoiso, int numAdducts, 
                              int tolerance, int numMass, NumericVector R_monoisitopeMassVector, 
                              NumericVector R_adductMassVector, List R_isotopes,
                              NumericVector R_isotopeListOrder, NumericVector R_massAxis)
{
  //Fill the class data structure with the information of the experiment
  AdductPairer::AdductDef	myAdductDef;  
  myAdductDef.numMass = numMass;
  myAdductDef.numAdducts = numAdducts;
  myAdductDef.numMonoiso = numMonoiso;
  myAdductDef.tolerance = tolerance;
  
  List result;
  AdductPairer myAdductPairer(&myAdductDef, R_monoisitopeMassVector, R_adductMassVector, R_isotopes, R_isotopeListOrder, R_massAxis); //Class constructor
  

  
  result = myAdductPairer.Run();  //Program executioning function
  
  return result;
}


