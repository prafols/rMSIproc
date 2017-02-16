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
#include <iostream>
#include <fstream>
#include <sstream>
#include "peakmatrixio.h"
using namespace Rcpp;

PeakMatrixIO::PeakMatrixIO(String dir_path):
path(dir_path)
{
  int nimg, ncol, nrow;
  int totalRows = 0;
  String defFile = path;
  defFile += "/data.def";
  std::ifstream file(defFile.get_cstring(), std::ios::in|std::ios::binary);
  if (file.is_open())
  {
    file.read ( (char*)(&ncol), sizeof(int));
    file.read ( (char*)(&nimg), sizeof(int));
    nRows = new IntegerVector(nimg);
    for(int i = 0; i < nimg; i++)
    {
      file.read ( (char*)(&nrow), sizeof(int));
      (*nRows)[i] = nrow;
      totalRows += nrow;
    }
    
    intMat = new NumericMatrix(totalRows, ncol);
    snrMat = new NumericMatrix(totalRows, ncol);
    areaMat = new NumericMatrix(totalRows, ncol);
    massVec = new NumericVector(ncol);
    posMat = new IntegerMatrix(totalRows, 2);
    sNames = new StringVector(nimg);
    
    //Reading names.txt
    String nameFile = path;
    nameFile += "/names.txt";
    std::ifstream fileNames(nameFile.get_cstring(), std::ios::in);
    std::string str; 
    if (fileNames.is_open())
    {
      int ip = 0;
      while (std::getline(fileNames, str))
      {
        // Process str that contains each line of names.txt 
        (*sNames)[ip] = str;
        ip++;
      }
      fileNames.close();
      if( ip != nimg)
      {
        stop("Fatal Error: number of samples names in names.txt does not match to number of images in peak matrix. Aborting...\n");
      }
    }
    else
    {
      Rcout <<"Warning: No names.txt file found creating unknown names...\n"; 
      for( int i = 0; i < nimg; i++)
      {
        (*sNames)[i] = "unknown sample";
      }
    }
  }
}

PeakMatrixIO::~PeakMatrixIO()
{
  delete intMat;
  delete snrMat;
  delete areaMat;
  delete massVec;
  delete posMat;
  delete nRows;
  delete sNames;
}

List PeakMatrixIO::LoadPeakMatrix()
{
  Rcout <<"Loading intensity matrix\n";
  LoadMat(intensity);
  Rcout <<"Loading SNR matrix\n";
  LoadMat(SNR);
  Rcout <<"Loading area matrix\n";
  LoadMat(area);
  Rcout <<"Loading mass vector\n";
  LoadMass();
  Rcout << "Loading pixel positions\n";
  LoadPos();
  return List::create( Named("mass") = *massVec, Named("intensity") = *intMat, Named("SNR") = *snrMat, Named("area") = *areaMat,
                             Named("pos") = *posMat, Named("numPixels") = *nRows, Named("names") = *sNames);
}

void PeakMatrixIO::StorePeakMatrix(List lpeak)
{
  int nrow = as<NumericMatrix>(lpeak["intensity"]).rows();
  int ncol = as<NumericMatrix>(lpeak["intensity"]).cols();
  int nImg = as<IntegerVector>(lpeak["numPixels"]).length();
  
  intMat = new NumericMatrix(nrow, ncol);
  snrMat = new NumericMatrix(nrow, ncol);
  areaMat = new NumericMatrix(nrow, ncol);
  massVec = new NumericVector(ncol);
  posMat = new IntegerMatrix(nrow, 2);
  nRows = new IntegerVector(nImg);
  sNames = new StringVector(nImg);
  
  *intMat = as<NumericMatrix>(lpeak["intensity"]);
  *snrMat = as<NumericMatrix>(lpeak["SNR"]);
  *areaMat = as<NumericMatrix>(lpeak["area"]);
  *massVec = as<NumericVector>(lpeak["mass"]);
  *posMat = as<IntegerMatrix>(lpeak["pos"]);
  *nRows = as<IntegerVector>(lpeak["numPixels"]);
  *sNames = as<StringVector>(lpeak["names"]);
  
  Rcout <<"Storing intensity matrix\n";
  StoreMat(intensity);
  Rcout <<"Storing SNR matrix\n";
  StoreMat(SNR);
  Rcout <<"Storing area matrix\n";
  StoreMat(area);
  Rcout <<"Storing mass vector\n";
  StoreMass();
  Rcout <<"Storing pixel positions\n";
  StorePos();
  
  Rcout <<"Storing definition file\n";
  String defFile = path;
  defFile += "/data.def";
  std::ofstream file(defFile.get_cstring(), std::ios::out|std::ios::binary|std::ios::trunc);
  if (file.is_open())
  {
    file.write((char*)(&ncol), sizeof(int)); 
    file.write((char*)(&nImg), sizeof(int));
    for(int i = 0; i < nImg; i++)
    {
      nrow = (*nRows)[i];
      file.write((char*)(&nrow), sizeof(int));
      
    }
    file.close();
  }
  else
  {
    stop("C Writing error: Def. file can not be opened.\n");
  } 
  
  //Writing the names.txt file
  String namesFile = path;
  namesFile += "/names.txt";
  std::ofstream fileNames(namesFile.get_cstring(), std::ios::out|std::ios::trunc);
  if(fileNames.is_open())
  {
    for(int i = 0; i < nImg; i++)
    {
      fileNames <<  ((*sNames)[i]) << "\n";
    }
    fileNames.close();
  }
  else
  {
    stop("C Writing error: Names file can not be opened.\n");
  } 
  
}

void PeakMatrixIO::StoreMat(DataType mt)
{
  NumericMatrix *mat = setMatPointer(mt);
  double buffer[mat->cols()]; //Writing  buffer
  std::ofstream file(getFileName(mt).get_cstring(), std::ios::out|std::ios::binary|std::ios::trunc);
  if (file.is_open())
  {
    for( int i = 0; i < mat->rows(); i++)
    {
      //Copy data to buffer
      for( int j = 0; j  < mat->cols(); j++)
      {
        buffer[j] = (*mat)(i, j);
      }
      file.write((char*)buffer, sizeof(double)*mat->cols()); 
    }
    file.close();
  }
  else
  {
    stop("C Writing error: File can not be opened.\n");
  }
}

void PeakMatrixIO::LoadMat(DataType mt)
{
  NumericMatrix *mat = setMatPointer(mt);
  char buffer[sizeof(double)*mat->cols()]; //Reading  buffer
  std::ifstream file(getFileName(mt).get_cstring(), std::ios::in|std::ios::binary);
  if (file.is_open())
  {
    for( int i = 0; i < mat->rows(); i++)
    {
      file.read (buffer, sizeof(double)*mat->cols());  
      
      //Copy data to mat
      for( int j = 0; j  < mat->cols(); j++)
      {
        (*mat)(i, j) = *((double*)( buffer + j*sizeof(double) ));
      }
    }
    file.close();
  }
  else
  {
    stop("C Reading error: File can not be open.\n");
  }
}

void PeakMatrixIO::StoreMass()
{
  double buffer[massVec->length()]; //Writing  buffer
  std::ofstream file(getFileName(mass).get_cstring(), std::ios::out|std::ios::binary|std::ios::trunc);
  if (file.is_open())
  {
    //Copy data to buffer
    memcpy(buffer, massVec->begin(), sizeof(double)*massVec->length());
    file.write((char*)buffer, sizeof(double)*massVec->length()); 
    file.close();
  }
  else
  {
    stop("C Writing error: File can not be opened.\n");
  } 
}

void PeakMatrixIO::LoadMass()
{
  char buffer[sizeof(double)*massVec->length()]; //Reading  buffer
  std::ifstream file(getFileName(mass).get_cstring(), std::ios::in|std::ios::binary);
  if (file.is_open())
  {
    file.read (buffer, sizeof(double)*massVec->length());  
    memcpy(massVec->begin(), buffer, sizeof(double)*massVec->length());
    file.close();
  }
  else
  {
    stop("C Reading error: File can not be open.\n");
  }
}

void PeakMatrixIO::StorePos()
{
  int buffer[posMat->cols()]; //Writing  buffer
  std::ofstream file(getFileName(pos).get_cstring(), std::ios::out|std::ios::binary|std::ios::trunc);
  if (file.is_open())
  {
    for(int i = 0; i < posMat->rows(); i ++)
    {
      for(int j = 0; j < posMat->cols(); j++)
      {
        buffer[j] = (*posMat)(i, j);
      }
      file.write((char*)buffer, sizeof(int)*posMat->cols()); 
    }
    file.close();
  }
  else
  {
    stop("C Writing error: File can not be opened.\n");
  } 
}

void PeakMatrixIO::LoadPos()
{
  char buffer[sizeof(int)*posMat->cols()]; //Reading  buffer
  std::ifstream file(getFileName(pos).get_cstring(), std::ios::in|std::ios::binary);
  if (file.is_open())
  {
    for(int i = 0; i < posMat->rows(); i ++)
    {
      file.read(buffer, sizeof(int)*posMat->cols());  
      for(int j = 0; j < posMat->cols(); j++)
      {
        (*posMat)(i, j) = *((int*)( buffer + j*sizeof(int) ));;
      }
    }
    file.close();
  }
  else
  {
    stop("C Reading error: File can not be open.\n");
  } 
}

String PeakMatrixIO::getFileName(DataType mt)
{
  //TODO: improve that method to avoid / char as path separator
  String fname = path;
  switch(mt)
  {
  case intensity:
    fname += "/peakint.mat";
    break;
  case SNR:
    fname += "/peaksnr.mat";
    break;
  case area:
    fname += "/peakarea.mat";
    break;
  case mass:
    fname += "/peakmass.vec";
    break;
  case pos:
    fname += "/pixelpos.mat";
    break;
  default:
    stop("Error: Invalid MatType mode.\n");
  }
  return fname;
}

NumericMatrix *PeakMatrixIO::setMatPointer(DataType mt)
{
  NumericMatrix *mat;
  switch(mt)
  {
  case intensity:
    mat = this->intMat;
    break;
  case SNR:
    mat = this->snrMat;
    break;
  case area:
    mat = this->areaMat;
    break;
  default:
    stop("Error: invalid DataType.\n");
  }
  return mat;
}

//// R exports /////////////////////////////////////////
//'LoadPeakMatrix.
//'
//'Loads a binned peaks matrix from HDD.
//'
//'@param path full path to directory from where data must be loaded.
//'@return  an R List containing intensity, SNR and area matrices and mass axis vector.
//'
// [[Rcpp::export]]
List LoadPeakMatrixC( String path )
{
  PeakMatrixIO pkMatObj(path);
  return pkMatObj.LoadPeakMatrix();
}

//'StorePeakMatrix.
//'
//'Stores a binned peaks matrix to HDD.
//'
//'@param path full path to directory where data must be stored.
//'@param mat an R List containing intensity, SNR and area matrices and mass axis vector.
//'
// [[Rcpp::export]]
void StorePeakMatrixC( String path, List mat )
{
  PeakMatrixIO pkMatObj(path);
  pkMatObj.StorePeakMatrix(mat);
}
