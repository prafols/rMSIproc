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


AdductPairer::AdductPairer(AdductDef *imgRunInfo, NumericVector R_monoisitopeMassVector, NumericVector R_adductMassVector)
{
	RunDef = imgRunInfo;

	monoisitopeMassVector = new double[RunDef->numMonoiso];
	for(int i = 0; i < RunDef->numMonoiso; i++)
	{
		monoisitopeMassVector[i] = R_monoisitopeMassVector[i];
	}

	mdVlength = RunDef->numMonoiso-1;
	maxShift = RunDef->numMonoiso-1;
	massDifferencesMatrix = new double*[maxShift];
  for(int i = 0; i < maxShift; i++)
  {
    massDifferencesMatrix[i] = new double[maxShift];
    for(int j = 0; j < maxShift; j++)
    {
      massDifferencesMatrix[i][j] = 0;
    }
  }

  adductMaskMatrix = new int*[maxShift];
  for(int i = 0; i < maxShift; i++)
  {
    adductMaskMatrix[i] = new int[maxShift];
    for(int j = 0; j < maxShift; j++)
    {
      adductMaskMatrix[i][j] = -1;	//-1 means that the combination is not related with adducts of any type
    }
  }

	adductMassVector = new double[RunDef->numAdducts];
	for(int i = 0; i < RunDef->numAdducts; i++)
	{
		adductMassVector[i] = R_adductMassVector[i];
	}

	currentShift = 0;
}



AdductPairer::~AdductPairer() //TODO all delets pls!!
{
	delete[] monoisitopeMassVector;
  
  for(int i = 0; i < maxShift; i++)
  {
	  delete[] massDifferencesMatrix[i];
  }
  delete[] massDifferencesMatrix;
  
  for(int i = 0; i < maxShift; i++)
  {
    delete[] adductMaskMatrix[i];
  }
  delete[] adductMaskMatrix;
  
	delete[] adductMassVector;
	
	delete[] adductMassDifferencesVector;
}



void AdductPairer::ShiftAndSubstract()
{
	currentShift++; 
	mdVlength = RunDef->numMonoiso-currentShift;
	for(int i = 0; i < mdVlength; i++)
	{
		massDifferencesMatrix[currentShift-1][i] = monoisitopeMassVector[currentShift+i] - monoisitopeMassVector[i];
	}
}



void AdductPairer::CalculateAdductMassDifferences()
{
	adductDiffCombinations = 0;	
	for(int i = 0; i < RunDef->numAdducts; i++)
	{
		adductDiffCombinations += i;
	}
	adductMassDifferencesVector = new double[adductDiffCombinations];

	int cnt1 = 0;
	int cnt2 = 1;
	int lim = RunDef->numAdducts-1;

	for(int i = 0; i < adductDiffCombinations; i++)
	{
		adductMassDifferencesVector[i] = adductMassVector[cnt1]-adductMassVector[cnt2];
		cnt2++;
		if(cnt2 >lim)
		{
			cnt1++;
			cnt2 = cnt1+1;
		}
	}
}



void AdductPairer::CheckDifferences()
{
	double ppmMasError = 0; //ppm mass error for each pair

	for(int i = 0; i < adductDiffCombinations; i++)
	{
		for(int j = 0; i < mdVlength; j++)
		{
			ppmMasError = (1000000*abs(massDifferencesMatrix[currentShift-1][j]-adductMassDifferencesVector[i]))/adductMassDifferencesVector[i];
			if(ppmMasError < RunDef->tolerance)
			{
				adductMaskMatrix[currentShift-1][j] = i; //Save the adduct combination of the candidate
				massDifferencesMatrix[currentShift-1][j] = ppmMasError;	//Save the ppm error
			}	else 
				{
					massDifferencesMatrix[currentShift-1][j] = -1;
				}
		}
	}
}

List AdductPairer::GenerateResultsList()
{
	List Results(2);
	NumericMatrix mdM(maxShift,maxShift);
	NumericMatrix amm(maxShift,maxShift);

	for(int i = 0; i< maxShift; i++)
	{
		for(int j = 0; j< maxShift; j++)
		{
			mdM(i,j) = massDifferencesMatrix[i][j];
			amm(i,j) = adductMaskMatrix[i][j];
		}
	}

	Results[0] = mdM(Range(0,maxShift),Range(0,maxShift));
	Results[1] = amm(Range(0,maxShift),Range(0,maxShift));

}


List AdductPairer::Run()
{
	CalculateAdductMassDifferences();

	while(currentShift<maxShift)
	{
		ShiftAndSubstract();
		CheckDifferences();
	}

	List Results = GenerateResultsList();
	return(Results);
}


// [[Rcpp::export]]
Rcpp::List adductAnnotation(int numMonoiso, int numAdducts, int tolerance, 
														NumericVector R_monoisitopeMassVector, NumericVector R_adductMassVector)
{
  //Fill the class data structure with the information of the experiment
  AdductPairer::AdductDef	myAdductDef;  
  myAdductDef.numAdducts = numAdducts;
  myAdductDef.numMonoiso = numMonoiso;
  myAdductDef.tolerance = tolerance;
  
  List result;
  AdductPairer myAdductPairer(&myAdductDef, R_monoisitopeMassVector, R_adductMassVector); //Class constructor
  result = myAdductPairer.Run();  //Algorithm run
  
  return result;
}


