/*************************************************************************
 *     rMSIproc - R package for MSI data processing
 *     Copyright (C) 2014 Pere Rafols Soler
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

#ifndef PEAK_MATRIX_IO_H
  #define PEAK_MATRIX_IO_H

#include <Rcpp.h>

class PeakMatrixIO
{
  public:
    PeakMatrixIO(Rcpp::String dir_path);
    ~PeakMatrixIO();
    
    Rcpp::List LoadPeakMatrix();
    void StorePeakMatrix(Rcpp::List lpeak);

  private:
    Rcpp::String path;
    enum DataType {intensity, SNR, area, mass, pos, motors};
    Rcpp::NumericMatrix *intMat;
    Rcpp::NumericMatrix *snrMat;
    Rcpp::NumericMatrix *areaMat;
    Rcpp::NumericVector *massVec;
    Rcpp::IntegerMatrix *posMat;
    Rcpp::IntegerMatrix *motorMat;
    Rcpp::IntegerVector *nRows;
    Rcpp::StringVector *sNames;
  
    void StoreMat(DataType mt);
    void LoadMat(DataType mt);
    void StoreMass();
    void LoadMass();
    void StorePos(DataType mt);
    void LoadPos(DataType mt);
    Rcpp::String getFileName(DataType mt);
    Rcpp::NumericMatrix *setMatPointer(DataType mt);
    Rcpp::IntegerMatrix *setPosPointer(DataType mt);
    Rcpp::DataFrame LoadNormalizationMatrix();
    void StoreNormalizationMatrix(Rcpp::DataFrame norms);
};

#endif