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
    int numPixels;
  }AdductDef;
  
  AdductPairer(AdductDef *imgRunInfo, NumericVector R_monoisitopeMassVector,
               NumericVector R_adductMassVector, List R_isotopes,
               NumericVector R_isotopeListOrder, NumericVector R_massAxis,
               NumericMatrix R_peakMatrix);  //Constructor
  
  ~AdductPairer();      //Destructor
  List Run();           //Program executioning function

private:
	AdductDef *RunDef;                    //Pointer to the structure.
	double *monoisitopeMassVector;        //Vector containing the masses of the monoisotopic species (length = numMonoiso)
	double **ppmMatrix;                   //Vector containing the substraction of the monoisotopic masses. At each step, the length is reduced by 1. (length = mdVlength)  
  double **massMatrix;
  double **peakMatrix;
  double *adductMassVector; 		        //Vector containing the adduct forming chemical elements masses.   
  double *adductPairFirstMassVector;    //Vector containing the mass of the first adduct of the pair. Ex. (K & Na , K & H, Na & H) -> (K, K, Na) 
  double *adductPairSecondMassVector;   //Vector containing the mass of the first adduct of the pair. Ex. (K & Na , K & H, Na & H) -> (Na, H, H) 
  double *massAxis;                     //Vector containint the mass axis of the peak matrix.
  int **adductPairsNameMatrix;          //Matrix containing the adduct elements of the pairs. -1 if no pair found.
  int *isotopeListOrder;                //Vector containing where the ion is found in the isotopes list
  int currentShift;							        //Current shift in the mass substraction
  int maxShift;                         //Maximum allowed shift for the number of monoisotopic masses
  int adductDiffCombinations;           //Length of adductMassDifferencesVector
  int mdVlength; 							          //Length of the mass differences vector. At each shift
  int positiveTest;                     //Number of pairs that have succes in the test
  List isotopes;                        //Results from the isotope test
  
  void AdductPairCreator();                   //Fills adductPairFirstMassVector and adductPairSecondMassVector
	void ShiftandCalculateAdductPairs();        //Substracts the monoisotopeMassVector with the shifted version of himself
  void ValidateCandidates();                  //Reads the information form the isotopes tests and merge it with the adduct candidates information.
  List GenerateResultsList();                 //Output construction function
};



double correlation(NumericVector v1, NumericVector v2, int vector_length)
{
  double sumXY = 0, sumX = 0, sumY = 0, sumX2 = 0, sumY2 = 0, num = 0, den = 0, result = 0;
  
  for(int i = 0; i < vector_length ; i++)
  {
    sumXY  += v1[i]*v2[i];
    sumX   += v1[i];
    sumY   += v2[i];
    sumX2  += v1[i]*v1[i];
    sumY2  += v2[i]*v2[i];
  }
  den = (vector_length*(sumX2)-(sumX*sumX))*(vector_length*(sumY2)-(sumY*sumY));
  den = sqrt(den);
  num = (vector_length*sumXY-(sumX*sumY));
  result = num/den;
  return(result);
}












