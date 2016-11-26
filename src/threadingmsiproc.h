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

#ifndef MULTI_THREAD_WORKER_H
  #define MULTI_THREAD_WORKER_H

#include <Rcpp.h>
#include <boost/thread.hpp>
#include "rmsicdataio.h"

class ThreadingMsiProc
{
  public:
    ThreadingMsiProc(int numberOfThreads, bool overWriteRamdisk, Rcpp::String basePath, Rcpp::StringVector fileNames, int massChannels, int *numRows, Rcpp::String dataType );
    ~ThreadingMsiProc();
    
  protected:
    //Pure virtual function to be implemented in ThreadingMsiProc class derivations.
    virtual void ProcessingFunction(int threadSlot);
    
    //Function to control threaded execution
    void runMSIProcessingCpp();
    
    CrMSIDataIO::DataCube **cubes; //Array of data cubes pointer, the length of this array will be the number of processing threads.
    
  private:  
    //The function to be executed for each thread. 
    //Data will be accessed from each thread using in-class member data and the index provided as threadSlot parameter.
    void ProcessingThread( int threadSlot );
    
    //Function to suspend main thread execution until some thread ends.
    void WaitForSomeThreadEnd();
    
    boost::mutex mtx; //Lock mechanism for signalling bDataReady vector
    int *iCube; //This vector will porvide which cube ID is processing each thread
    bool *bDataReady; //This vector will contain true when a worker thread completes a datacube processing
    
    //Condition variable to notify thread ends
    boost::condition_variable  life_end_cond;
    boost::mutex               life_end_mutex;
    bool                       life_end;
    int numOfThreads;
    CrMSIDataIO *ioObj; //Data access object
    boost::thread *tworkers; //Thread objects
    bool bDataOverWrite;
};
  
#endif