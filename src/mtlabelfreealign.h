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

#ifndef MT_ALIGN_H
  #define MT_ALIGN_H
#include <Rcpp.h>
#include <boost/thread.hpp>
#include "labelfreealign.h"
#include "threadingmsiproc.h"


class MTLabelFreeAlign : public ThreadingMsiProc 
{
  public:
    //Data structur used to completly define the processing pipeline
    typedef struct
    {
      double *mass; //pointer to the global mass axis
      double *ref_spectrum; //Reference spectrum used for alignment
      int massChannels; //Number of points in the mass axis array
      bool bilinearMode; //If the alignment must be performed linear or bilinear
      Rcpp::StringVector fileNames; //Filname of each ramdisk file
      int *numRows; //An array containing the number of rows stored in each ramdisk file. The length ot this array is the length of fileNames
      Rcpp::String dataType; //A string with the data type
      int AlignmentIterations;
      int AlignmentMaxShift;
      double RefLow;
      double RefMid;
      double RefHigh;
      int OverSampling;
      int numOfThreads;
    }ImgProcDef;  
    
    MTLabelFreeAlign(ImgProcDef imgRunInfo);
    ~MTLabelFreeAlign();
    
    //Exectur a full imatge processing using threaded methods and returns the used shifts in first iteration
    Rcpp::List Run(); 

  private:
    LabelFreeAlign **alngObj;
    int numOfPixels;
    LabelFreeAlign::TLags *mLags; //A place to store alignment lags

    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
    
    boost::mutex fftSharedMutex;
};
#endif
