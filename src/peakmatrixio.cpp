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
    motorMat = new IntegerMatrix(totalRows, 2);
    sNames = new StringVector(nimg);
    sUUID = new StringVector(nimg);
    
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
    
    //Reading uuid.txt
    String uuidFile = path;
    uuidFile += "/uuid.txt";
    std::ifstream fileUUIDs(uuidFile.get_cstring(), std::ios::in);
    if (fileUUIDs.is_open())
    {
      int ip = 0;
      while (std::getline(fileUUIDs, str))
      {
        // Process str that contains each line of names.txt 
        (*sUUID)[ip] = str;
        ip++;
      }
      fileUUIDs.close();
      if( ip != nimg)
      {
        stop("Fatal Error: number of samples names in uuid.txt does not match to number of images in peak matrix. Aborting...\n");
      }
    }
    else
    {
      Rcout <<"Warning: No uuid.txt file found creating empty UUID's...\n"; 
      for( int i = 0; i < nimg; i++)
      {
        (*sUUID)[i] = "";
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
  delete motorMat;
  delete nRows;
  delete sNames;
  delete sUUID;
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
  LoadPos(pos);
  
  
  List retList = List::create( Named("mass") = *massVec, Named("intensity") = *intMat, Named("SNR") = *snrMat, Named("area") = *areaMat,
                                 Named("pos") = *posMat, Named("numPixels") = *nRows, Named("names") = *sNames, Named("uuid") = *sUUID );
  
  // check if motorpos.mat and load motor coords if so
  String motorTxtFile = path;
  motorTxtFile += "/motorpos.mat";
  std::ifstream fileMotorsTxt(motorTxtFile.get_cstring());
  if (fileMotorsTxt.is_open())
  {
    Rcout << "Loading motors coordinates matrix\n";
    LoadPos(motors);
    retList["posMotors"] = *motorMat;
    fileMotorsTxt.close();
  }
  
  // check if norms.txt and load normalizations if so
  String normsTxtFile = path;
  normsTxtFile += "/norms.txt";
  std::ifstream fileNormsTxt(normsTxtFile.get_cstring());
  if (fileNormsTxt.is_open())
  {
    Rcout << "Loading normalizations\n";
    retList["normalizations"] = LoadNormalizationMatrix();
    fileNormsTxt.close();
  }
  
  return retList;
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
  motorMat = new IntegerMatrix(nrow, 2);
  nRows = new IntegerVector(nImg);
  sNames = new StringVector(nImg);
  sUUID = new StringVector(nImg);
  
  *intMat = as<NumericMatrix>(lpeak["intensity"]);
  *snrMat = as<NumericMatrix>(lpeak["SNR"]);
  *areaMat = as<NumericMatrix>(lpeak["area"]);
  *massVec = as<NumericVector>(lpeak["mass"]);
  *posMat = as<IntegerMatrix>(lpeak["pos"]);
  if( lpeak.containsElementNamed("posMotors"))
  {
    *motorMat = as<IntegerMatrix>(lpeak["posMotors"]);
  }
  else
  {
    *motorMat = as<IntegerMatrix>(lpeak["pos"]);
  }
  *nRows = as<IntegerVector>(lpeak["numPixels"]);
  *sNames = as<StringVector>(lpeak["names"]);
  if( lpeak.containsElementNamed("uuid"))
  {
    *sUUID = as<StringVector>(lpeak["uuid"]);
  }
  else
  {
    *sUUID = StringVector(sNames->length()); //No found ID fields so init with empty uuid's
  }
  
  Rcout <<"Storing intensity matrix\n";
  StoreMat(intensity);
  Rcout <<"Storing SNR matrix\n";
  StoreMat(SNR);
  Rcout <<"Storing area matrix\n";
  StoreMat(area);
  Rcout <<"Storing mass vector\n";
  StoreMass();
  Rcout <<"Storing pixel positions\n";
  StorePos(pos);
  Rcout <<"Storing motor coordinates\n";
  StorePos(motors);
  
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
  
  //Writing the uuid.txt file
  String uuidFile = path;
  uuidFile += "/uuid.txt";
  std::ofstream fileUUIDs(uuidFile.get_cstring(), std::ios::out|std::ios::trunc);
  if(fileUUIDs.is_open())
  {
    for(int i = 0; i < nImg; i++)
    {
      fileUUIDs <<  ((*sUUID)[i]) << "\n";
    }
    fileUUIDs.close();
  }
  else
  {
    stop("C Writing error: UUID's file can not be opened.\n");
  } 
  
  //Check if there ara also normalizations to store
  if( lpeak.containsElementNamed("normalizations") )
  {
    Rcout <<"Storing normalizations files\n";
    DataFrame norm = as<DataFrame>(lpeak["normalizations"]);
    StoreNormalizationMatrix(norm);
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

void PeakMatrixIO::StorePos(DataType mt)
{
  IntegerMatrix *ptrMat = setPosPointer(mt);
  int buffer[ptrMat->cols()]; //Writing  buffer
  std::ofstream file(getFileName(mt).get_cstring(), std::ios::out|std::ios::binary|std::ios::trunc);
  if (file.is_open())
  {
    for(int i = 0; i < ptrMat->rows(); i ++)
    {
      for(int j = 0; j < ptrMat->cols(); j++)
      {
        buffer[j] = (*ptrMat)(i, j);
      }
      file.write((char*)buffer, sizeof(int)*ptrMat->cols()); 
    }
    file.close();
  }
  else
  {
    stop("C Writing error: File can not be opened.\n");
  } 
}

void PeakMatrixIO::LoadPos(DataType mt)
{
  IntegerMatrix *ptrMat = setPosPointer(mt);
  char buffer[sizeof(int)*ptrMat->cols()]; //Reading  buffer
  std::ifstream file(getFileName(mt).get_cstring(), std::ios::in|std::ios::binary);
  if (file.is_open())
  {
    for(int i = 0; i < ptrMat->rows(); i ++)
    {
      file.read(buffer, sizeof(int)*ptrMat->cols());  
      for(int j = 0; j < ptrMat->cols(); j++)
      {
        (*ptrMat)(i, j) = *((int*)( buffer + j*sizeof(int) ));;
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
  case motors:
    fname += "/motorpos.mat";
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

IntegerMatrix *PeakMatrixIO::setPosPointer(DataType mt)
{
  IntegerMatrix *mat;
  switch(mt)
  {
  case pos:
    mat = this->posMat;
    break;
  case motors:
    mat = this->motorMat;
    break;
  default:
    stop("Error: invalid DataType.\n");
  }
  return mat;
}

DataFrame PeakMatrixIO::LoadNormalizationMatrix()
{
  //1 Read the norms.txt file to get number of columns and normalization names.
  String normsTxtFile = path;
  normsTxtFile += "/norms.txt";
  std::ifstream fileNormNames(normsTxtFile.get_cstring(), std::ios::in);
  std::string str; 
  CharacterVector normNames;
  if (fileNormNames.is_open())
  {
    while (std::getline(fileNormNames, str))
    {
      normNames.push_back(str);
    }
    fileNormNames.close();
  }
  else
  {
    stop("Error: No norms.txt file found but this function was called, this is a bug...\n"); 
  }
  
  Rcout << "Found " << normNames.length() << " normalizations in norms.txt\n";
  
  List tmpNorms( normNames.length() );
  
  char buffer[sizeof(double)*posMat->rows()]; //Reading  buffer
  String normsMatFile = path;
  normsMatFile += "/norms.mat";
  std::ifstream file(normsMatFile.get_cstring(), std::ios::in|std::ios::binary);
  if (file.is_open())
  {
    for( int i = 0; i < tmpNorms.length(); i++)
    {
      file.read (buffer, sizeof(double)*posMat->rows());  
      tmpNorms[i] = NumericVector(posMat->rows());
      memcpy(as<NumericVector>(tmpNorms[i]).begin(), buffer, sizeof(double)*(posMat->rows()));
    }
    file.close();
  }
  else
  {
    stop("C Reading error: File can not be open.\n");
  }
  
  //Create  the dataframe from the tmp list and names
  DataFrame DFnorms(tmpNorms);
  DFnorms.attr("names") = normNames;

  return(DFnorms);
}

void PeakMatrixIO::StoreNormalizationMatrix(DataFrame norms)
{
  //Store normalization data.frame streaming data by columns (each column is a noramalization vector)
  double buffer[norms.nrows()]; //Writing  buffer
  NumericVector dfColumn;
  String fname = path;
  fname += "/norms.mat";
  std::ofstream file(fname.get_cstring(), std::ios::out|std::ios::binary|std::ios::trunc);
  if (file.is_open())
  {
    for( int i = 0; i < norms.size(); i++)
    {
      //Copy data to buffer
      dfColumn = norms[i];
      memcpy(buffer, dfColumn.begin(), sizeof(double)*dfColumn.length());
      file.write((char*)buffer, sizeof(double)*norms.nrows()); 
    }
    file.close();
  }
  else
  {
    stop("C Writing error: File can not be opened.\n");
  }
  
  //Store the names in a plain text file
  String namesFile = path;
  namesFile += "/norms.txt";
  std::ofstream fileNames(namesFile.get_cstring(), std::ios::out|std::ios::trunc);
  StringVector normNames = norms.names();
  if(fileNames.is_open())
  {
    for(int i = 0; i < normNames.length(); i++)
    {
      fileNames <<  normNames[i] << "\n";
    }
    fileNames.close();
  }
  else
  {
    stop("C Writing error: Names file can not be opened.\n");
  } 
}

//// R exports /////////////////////////////////////////
//'LoadPeakMatrix.
//'
//'Loads a binned peaks matrix from HDD.
//'
//'@param path full path to directory from where data must be loaded.
//'@return  an R List containing intensity, SNR and area matrices, mass axis vector and if available the normalizations data.frame.
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
//'@param mat an R List containing intensity, SNR and area matrices the mass axis vector and an R data.frame containing a normalization on each column.
//'
// [[Rcpp::export]]
void StorePeakMatrixC( String path, List mat )
{
  PeakMatrixIO pkMatObj(path);
  pkMatObj.StorePeakMatrix(mat);
}
