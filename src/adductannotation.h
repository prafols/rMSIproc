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
#include <Rcpp.h>
#include <stdlib.h>
#include <cmath>

using namespace Rcpp;

class AdductPairer
{
public:
  typedef struct 
  {
    int numMonoiso; //Number of monoistopic peaks.
    int numAdducts;	//Number of adduct forming chemical masses.  
    int tolerance;  //Mass tolerance in ppm
  }AdductDef;
  
  AdductPairer(AdductDef *imgRunInfo, NumericVector R_monoisitopeMassVector, NumericVector R_adductMassVector);  //Constructor
  ~AdductPairer();                      //Destructor
  List Run();                           //Program executioning function

private:
	AdductDef *RunDef;
	int mdVlength; 							          //Length of the mass differences vector. At each shift
	double *monoisitopeMassVector;          //Vector containing the masses of the monoisotopic species (length = numMonoiso)
	double **massDifferencesMatrix;         //Vector containing the substraction of the monoisotopic masses. At each step, the length is reduced by 1. (length = mDV.length)  
  double *adductMassVector; 		          //Vector containing the adduct forming chemical elements masses.   
  double *adductMassDifferencesVector;    //Vector containing the mass differeneces between adduct forming chemical elements.
  int **adductMaskMatrix;                //Boolean matrix with the same size of massDifferencesMatrix. TRUE means that the a pair of ians could be adducts. At the end, TRUE positions will write in the massDifferencesMatrix the ppm error and in the FALSE a -1 
  int currentShift;							          //Current shift in the mass substraction
  int maxShift;                           //Maximum allowed shift for the number of monoisotopic masses
  int adductDiffCombinations;             //Length of adductMassDifferencesVector
  
  void CalculateAdductMassDifferences();  //Fills adductMassDifferencesVector
	void ShiftAndSubstract();               //Substracts the monoisotopeMassVector with the shifted version of himself
  void CheckDifferences();                //Checks if there's any pair of peaks with the current
  List GenerateResultsList();
};