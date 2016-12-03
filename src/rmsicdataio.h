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
#ifndef RMSI_DATA_IO_H
  #define RMSI_DATA_IO_H
#include <Rcpp.h>

class CrMSIDataIO
{
  public: 
    //basePath: full path to the directory containing ramdisk files in *.dat format.
    //fileNames: a vector of names of each ramdisk files in basepath.
    //massChannels: the length of mass axis so the number of columns in the datacube.
    //numRows: a vector of integers with the same length as fileNames containing the number of rows of each datacube.
    //dataType: a string specifing the C data type used to store data in the ramdisk.
    CrMSIDataIO( Rcpp::String basePath, Rcpp::StringVector fileNames, int massChannels, int *numRows, Rcpp::String dataType );
    ~CrMSIDataIO();
    void printDataInfo();
    
    //Struct to define a whole data cube in memory
    typedef struct
    {
      int cubeID;
      int ncols;
      int nrows;
      double **data;  
    } DataCube;
    
    //Loads a data cube specified by iCube into data_ptr
    //WARNING: This is not a thread-safe method, it must be only used on the main thread who'll take care of loading data from HDD and copy it to other thread mem space.
    //It returns a pointer to an allocated structure containing the datacube.
    //It is responisability of user to free memmory of the loaded datacube using freeDataCube() function.
    DataCube *loadDataCube( int iCube);
    void freeDataCube(DataCube *data_ptr);
    
    //Stores a datacube to the path assosiated with its ID
    void storeDataCube(int iCube, DataCube *data_ptr);
    
    //Retunrs the total number of cubes in the ramdisk
    int getNumberOfCubes();
    
    //Returns the id of the first spectrum in the given cube
    int getFirstSpectrumIdInCube(int cube_id);

  private:
    Rcpp::String dataPath;
    Rcpp::StringVector ffFiles;
    int dataLength;
    int *rowCounts;
    Rcpp::String ffDataType;
    
    //Returns the full path to a ramdisk
    Rcpp::String getFullPath(int cube_id);
};
  
#endif