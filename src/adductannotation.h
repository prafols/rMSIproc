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
    int numMass;    //Length of the mass axis.
    int numMonoiso; //Number of monoistopic peaks.
    int numAdducts;	//Number of adduct forming chemical masses.  
    int tolerance;  //Mass tolerance in ppm
  }AdductDef;
  
  AdductPairer(AdductDef *imgRunInfo, NumericVector R_monoisitopeMassVector, NumericVector R_adductMassVector, List R_isotopes, NumericVector R_isotopeListOrder, NumericVector R_massAxis);  //Constructor
  ~AdductPairer();                        //Destructor
  List Run();                             //Program executioning function

private:
	AdductDef *RunDef;                          //Pointer to the structure.
	double *monoisitopeMassVector;              //Vector containing the masses of the monoisotopic species (length = numMonoiso)
	double **ppmMatrix;             //Vector containing the substraction of the monoisotopic masses. At each step, the length is reduced by 1. (length = mdVlength)  
  double *adductMassVector; 		              //Vector containing the adduct forming chemical elements masses.   
  double *adductPairFirstMassVector;          //Vector containing the mass of the first adduct of the pair. Ex. (K & Na , K & H, Na & H) -> (K, K, Na) 
  double *adductPairSecondMassVector;          //Vector containing the mass of the first adduct of the pair. Ex. (K & Na , K & H, Na & H) -> (Na, H, H) 
  double *massAxis;                           //Vector containint the mass axis of the peak matrix.
  int **adductPairsNameMatrix;                //Matrix containing the adduct elements of the pairs. -1 if no pair found.
  int *isotopeListOrder;                      //Vector containing where the ion is found in the isotopes list
  int currentShift;							              //Current shift in the mass substraction
  int maxShift;                               //Maximum allowed shift for the number of monoisotopic masses
  int adductDiffCombinations;                 //Length of adductMassDifferencesVector
  int mdVlength; 							                //Length of the mass differences vector. At each shift
  int positiveTest;                           //Number of pairs that have succes in the test
  List isotopes;                              //Results from the isotope test
  
  void CalculateAdductMassDifferences();      //Fills adductMassDifferencesVector
	void ShiftandCalculateAdductPairs();        //Substracts the monoisotopeMassVector with the shifted version of himself
  void ValidateCandidates();                  //Reads the information form the isotopes tests and merge it with the adduct candidates information.
  List GenerateResultsList();                 //Output construction function
};



