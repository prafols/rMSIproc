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
    int massPeaks;        //Number of mass Peaks
    long int massChannels;     //Number of mass channels
    int numPixels;        //Number of pixels
    int NumIso;           //Number of isotopes to be found given a specific monoisotopic mass
    int scansT;           //scans tolerance for peak candidates
    int CrntNumIso;       //Current isotop beeing computed
    int MaxCan;           //maximum number of candidates for one peak
    double accu;          //Score needed to pass the isotope test
  }IsoDef;
  
  Deisotoper(IsoDef imgRunInfo, double** PeakMtx, double* ProMassAxis, double* RawMassAxis);
  ~Deisotoper();
  
  //Program executioning function
  List Run(double **PeakMatrix, IsoDef imgRuninfo, double* massAxis, double* RawMassAxis);
  
private:
  double *pTPI;         //pointer to the total pixel intensity vector
  int *pPeakOrder;      //vector containing the columns ordered by total pixel intensity
  int *pRawPeakOrder;   //vector containing the index of the ordered columns by TPI in the raw mass axis
  int *pNumCandidates;  //vector containing the number of peak candidates for each selected mass
  int maxCan;           ////maximum number of candidates for one peak
  
  
  //private functions
  void PeakOrganizer(double** PeakMtx, IsoDef imgRuninfo, double* ProMassAxis, double* RawMassAxis);   //Selects the peaks by intensity and points each preocessed mass to the raw mass axis
  NumericMatrix CandidateFinder(double* ProMassAxis, double* RawMassAxis, IsoDef imgRuninfo, NumericVector PeaksToCheck);   //Finds the isotope candidates for each selected peak
  double* ScoreCalculator(int* test, int NumCan, double** PeakMatrix, IsoDef imgRuninfo, double* result, double* ProMassAxis, double before_m_slope);    //Computes the morphology & intensity test
  List MatrixAnnotator(double **PeakMatrix, IsoDef imgRuninfo, NumericMatrix CanMatrix, double* ProMassAxis, NumericVector PeaksToCheck);    //Annotates the matrix
};
#endif

