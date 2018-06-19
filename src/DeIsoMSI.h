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
#include <RcppGSL.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>

using namespace Rcpp;

class Deisotyper
{
  public:
    //Data structur used to completly define the processing pipeline
    typedef struct
    {
      int massChannels;     //Number of mass channels
      int numPixels;        //Number of pixels
      double IntPcent;      //Percentage of the most intense peak that a mass requires in order to enter the process of deisotyping
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
    int *pAnnotation;     //Vector containing the results from the test
    int *pNumCandidates;  //vector containing the number of peak candidates for each selected mass
    int FPOlength;        //length of the filtered pPeakOrderVector
    int FPOsubstractor;   //length of the checked pPeakOrderVector
    int MaxCan;           //maximum number of candidates for one peak
    
    //Constants from the lineal model that computes the intensity score
    const double c = 0.5; //standard deviation of the lineal model 0.05991173
    const double c0 = -0.129;  //origin correction of the lineal model
    const double c1 = 0.0008822;   //slope of the lineal model
      
    
    //private function
    void PeakFilter(double** PeakMtx, IsoDef imgRuninfo);           //Selects the peaks that are gonna be proccesed
    NumericMatrix PeakFinder(double* massAxis, IsoDef imgRuninfo);  //Finds the isotope candidates for each selected peak
    NumericMatrix CandidateMatrixCheck(NumericMatrix CanMatrix);    //Checks if there is any peak that it is already a candidate for another peak
    double* MorphologyTest(int* test, int NumCan, double** PeakMatrix, int numPixels, double* result);  //Computes the morphology test returning the morphology score
    double IntensityTest(int* test, int Candidate, double** PeakMatrix, int numPixels, double* massAxis);   //Computes the Intensity test returning the slope score
    double SlopeScore(double slope, double mass);   //Calculets the Slope socre
    NumericMatrix MatrixAnnotator(double **PeakMatrix, IsoDef imgRuninfo, NumericMatrix CanMatrix, double* massAxis);   //Annotates the matrix
};
#endif
