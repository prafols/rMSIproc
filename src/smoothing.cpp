/*************************************************************************
 * This file is part of rMSIproc.
 *
 * rMSIproc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * rMSIproc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rMSIproc.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#include <Rcpp.h>
#include <vector>
#include "smoothing.h"
using namespace Rcpp;

Smoothing::Smoothing(int SavitzkyGolayKerSize )
{
  switch(SavitzkyGolayKerSize)
  {
  case 5:
    sgC.push_back(-3.0/35.0);
    sgC.push_back(12.0/35.0);
    sgC.push_back(17.0/35.0);
    sgC.push_back(12.0/35.0);
    sgC.push_back(-3.0/35.0);
    break;
  case 7:
    sgC.push_back(-2.0/21.0);
    sgC.push_back(3.0/21.0);
    sgC.push_back(6.0/21.0);
    sgC.push_back(7.0/21.0);
    sgC.push_back(6.0/21.0);
    sgC.push_back(3.0/21.0);
    sgC.push_back(-2.0/21.0);
    break;
  case 9:
    sgC.push_back(-21.0/231.0);
    sgC.push_back(14.0/231.0);
    sgC.push_back(39.0/231.0);
    sgC.push_back(54.0/231.0);
    sgC.push_back(59.0/231.0);
    sgC.push_back(54.0/231.0);
    sgC.push_back(39.0/231.0);
    sgC.push_back(14.0/231.0);
    sgC.push_back(-21.0/231.0);
    break;
  case 11:
    sgC.push_back(-36.0/429.0);
    sgC.push_back(9.0/429.0);
    sgC.push_back(44.0/429.0);
    sgC.push_back(69.0/429.0);
    sgC.push_back(84.0/429.0);
    sgC.push_back(89.0/429.0);
    sgC.push_back(84.0/429.0);
    sgC.push_back(69.0/429.0);
    sgC.push_back(44.0/429.0);
    sgC.push_back(9.0/429.0);
    sgC.push_back(-36.0/429.0);
    break;
  default:
    stop("Error, not valid SavitzkyGolay kernel size, valid values are: 5, 7, 9 and 11");
  }
  
}

Smoothing::~Smoothing()
{
  
}

NumericVector Smoothing::smoothSavitzkyGolay(NumericVector x)
{
  double xC[x.length()];
  std::memcpy(xC, x.begin(), sizeof(double)*x.length());
  smoothSavitzkyGolay(xC, x.length());
  NumericVector y(x.length());
  std::memcpy(y.begin(), xC, sizeof(double)*y.length());
  return y;
}

//Performs the SavitzkyGolay smoothing and overwites the original data with smoothed data
void Smoothing::smoothSavitzkyGolay(double *x, int length)
{
  //Convolution with SavitzkyGolay kernel
  double y[length];
  for( int i = 0; i < length; i++)
  {
    y[i] = 0; //Init a zero
  }
  for( int i = (sgC.size() - 1)/2; i < length - ((sgC.size() - 1)/2); i++)
  {
    for( int j = 0; j < sgC.size(); j++)
    {
      y[i] = y[i] + x[i + j -((sgC.size() - 1)/2) ] * sgC[j];
    }
  }
  
  //Overwrite input pointer
  std::memcpy(x, y, sizeof(double)*length);
}


//' Smoothing_SavitzkyGolay.
//' 
//' Computes the Savitzky-Golay smoothing of a vector x using a filter size of sgSize.
//' @param x the data vector to smooth.
//' @param sgSize valid values are: 5, 7, 9 and 11.
//' @return the smoothed data vector.
// [[Rcpp::export]]
NumericVector Smoothing_SavitzkyGolay(NumericVector x, int sgSize = 5)
{
  Smoothing smObj(sgSize);
  return smObj.smoothSavitzkyGolay(x);
}


