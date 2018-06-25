/*************************************************************************
 * 
 *     Copyright (C) 2018 Lluc Sementé Fernàndez
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

#ifndef DEISOTYPER_H
  #define DEISOTYPER_H

#include <Rcpp.h>
#include <stdlib.h>
#include <cmath>


using namespace Rcpp;

class Deisotyper
{
  public:
    //Data structur used to completly define the processing pipeline
    typedef struct
    {
      int massChannels;     //Number of mass channels
      int numPixels;        //Number of pixels
      int NumIso;           //Number of isotopes to be found given a specific monoisotopic mass
      int ppmT;             //ppm tolerance for peak candidates
    }IsoDef;
    
    Deisotyper(IsoDef imgRunInfo, double** PeakMtx);
    ~Deisotyper();

    //Program executioning function
    NumericMatrix Run(double **PeakMatrix, IsoDef imgRuninfo, double* massAxis);
  
  private:
    double *pTPI;         //pointer to the total pixel intensity vector
    int *pPeakOrder;      //vector containing the columns ordered by total pixel intensity
    int *pNumCandidates;  //vector containing the number of peak candidates for each selected mass
    int MaxCan;           //maximum number of candidates for one peak
    
    //Constants from the lineal model that computes the intensity score
    const double c = 0.5; 
    const double c0 = -0.129;      //interception term
    const double c1 = 0.0008822;   //slope of the lineal model
      
    //private functions
    void PeakFilter(double** PeakMtx, IsoDef imgRuninfo);           //Selects the peaks that are gonna be proccesed
    NumericMatrix CandidateFinder(double* massAxis, IsoDef imgRuninfo);  //Finds the isotope candidates for each selected peak
    double* ScoreCalculator(int* test, int NumCan, double** PeakMatrix, int numPixels, double* result, double* massAxis);  //Computes the morphology & intensity test
    double SlopeScore(double slope, double mass);   //Calculets the Slope socre
    NumericMatrix MatrixAnnotator(double **PeakMatrix, IsoDef imgRuninfo, NumericMatrix CanMatrix, double* massAxis);   //Annotates the matrix
};
#endif
