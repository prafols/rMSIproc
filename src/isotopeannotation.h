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

#ifndef DEISOTOPER_H
  #define DEISOTOPER_H

#include <Rcpp.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace Rcpp;

class Deisotoper
{
  public:
    //Data structure used to completly define the algorithm run
    typedef struct
    {
      int massPeaks;             //Number of mass Peaks (length of the matrix mass axis)
      int massChannels;          //Number of mass channels (length of the image mass axis)
      int numPixels;             //Number of pixels
      int z;                     //Charge of the pattern to be found
      int numIso;                //Number of isotopes to be found given a specific monoisotopic mass
      int tolerance;             //Mass tolerance for peak candidates
      double scoreThreshold;     //Final score needed to pass the isotope test
      bool ToleranceInScans;     //If true, tolerance is applied in scans, else is in ppm
    }IsoDef;
    
    Deisotoper(IsoDef *imgRunInfo, NumericMatrix PeakMtx, NumericVector MatrixAxis, NumericVector ImageAxis); //Constructor
    ~Deisotoper();
    List Run();   //Program executioning function
    
  private:
    IsoDef *imgRuninfo;       //pointer to an IsoDef struct
    double **pMatrix;         //pointer to the peak matrix
    double *pMatrixAxis;      //pointer to the peak matrix mass axis
    double *pImageAxis;       //pointer to the image mass channels(scans) axis
    double *pTPI;             //pointer to the total pixel intensity vector
    double *massErrorCurve;   //Error curve in order to compute the mass error in ppm
    int *pMatrixPeakOrder;    //vector containing the columns ordered by total pixel intensity
    int *pImagePeakOrder;     //vector containing the index of the ordered columns by TPI in the raw mass axis
    int *pNumCandidates;      //vector containing the number of peak candidates for each selected mass
    int maxCan;               //Maximum number of candidates for one peak in the whole matrix
    int CrntNumIso;           //Current isotopic stage beeing evaluated

    void SortPeaksByIntensity();  //Sorts the peaks by intensity 
    void MatrixToImageAxis(); //Points the peak matrix mass axis to the image mass axis in order to use the tolerance in scans
    NumericMatrix CandidateFinder(NumericVector PeaksToCheck);  //Finds the isotope candidates for each selected peak
    double* ScoreCalculator(int *CandidateRow, int NumCan, double *result, double lastModslope, int peakNumber); //Computes all the scores over a candidates row 
    List MatrixAnnotator(NumericMatrix CanMatrix,NumericVector PeaksToCheck); //Directs the ScoreCalculator function for each ion to its candidates
    double getToleranceFromCurve(int massIndex);  //get the base tolerance for the mass error score. Scan mode only.
    double* getCandidateLimits(int massIndex); //Finds the candidates mass limits for each peak in both scan and ppm mode
};

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

#endif

