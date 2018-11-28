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
#include <cmath>

using namespace Rcpp;

class Deisotoper
{
public:
  //Data structur used to completly define the processing pipeline
  typedef struct
  {
    int massPeaks;             //Number of mass Peaks
    long int massChannels;     //Number of mass channels
    int numPixels;             //Number of pixels
    int numIso;                //Number of isotopes to be found given a specific monoisotopic mass
    int tolerance;             //Tolerance for peak candidates
    double scoreThreshold;             //Score needed to pass the isotope test
    bool ToleranceInScans;     //If true, tolerance is applied in scans, else is in ppm
  }IsoDef;
  
  Deisotoper(IsoDef *imgRunInfo, NumericMatrix PeakMtx, NumericVector MatrixAxis, NumericVector ImageAxis);
  Deisotoper(IsoDef *imgRunInfo, NumericMatrix PeakMtx, NumericVector ImageAxis);
  ~Deisotoper();
  
  //Program executioning function
  List Run(IsoDef *imgRuninfo);
  
private:
  double **pMatrix;
  double *pMatrixAxis;
  double *pImageAxis;
  double *pTPI;             //pointer to the total pixel intensity vector
  int *pMatrixPeakOrder;    //vector containing the columns ordered by total pixel intensity
  int *pImagePeakOrder;     //vector containing the index of the ordered columns by TPI in the raw mass axis
  int *pNumCandidates;      //vector containing the number of peak candidates for each selected mass
  int maxCan;               //Maximum number of candidates for one peak in the whole matrix (Just to not oversize the candidates matrix)
  int CrntNumIso;            //Current isotop beeing computed
  int MaxCan;                //Maximum number of candidates for one peak
  //double *massErrorCurve;  //TODO!!!! Error curve in order to compute the mass ppm error!
  
  //private functions
  void SortPeaksByIntensity(IsoDef *imgRuninfo);   //Selects the peaks by intensity and points each processed mass to the raw mass axis
  void MatrixToImageAxis(IsoDef *imgRuninfo);      //Points one mass axis to the other in order to use the tolerance in scans
  NumericMatrix CandidateFinder(IsoDef *imgRuninfo, NumericVector PeaksToCheck);   //Finds the isotope candidates for each selected peak
  double* ScoreCalculator(int* test, int NumCan,IsoDef *imgRuninfo, double* result, double before_m_slope);    //Computes the morphology & intensity test
  List MatrixAnnotator(IsoDef *imgRuninfo, NumericMatrix CanMatrix,NumericVector PeaksToCheck);    //Annotates the matrix
};
#endif

