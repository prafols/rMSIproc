/*************************************************************************
 * rMSIproc - R package for MSI data processing
 * Copyright (C) 2014 Pere Rafols Soler
 *
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
#ifndef SMOOTHING_H
  #define SMOOTHING_H

#include <Rcpp.h>
  
class Smoothing
{
  public:
    Smoothing(int SavitzkyGolayKerSize = 5);
    ~Smoothing();
    Rcpp::NumericVector smoothSavitzkyGolay(Rcpp::NumericVector x);
    void smoothSavitzkyGolay(double *x, int length);
    
  private:
    std::vector<double> sgC; //This is the SavitzkyGolay kernel
};


#endif