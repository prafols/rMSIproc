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

#include <Rcpp.h>
#include <boost/thread.hpp>
//#include <boost/bind.hpp>
#include "threadingmsiproc.h" 

ThreadingMsiProc::ThreadingMsiProc( int numberOfThreads, bool overWriteRamdisk, Rcpp::String basePath, Rcpp::StringVector fileNames, int massChannels, int *numRows, Rcpp::String dataType )
{
  numOfThreadsDouble = 2*numberOfThreads;
  ioObj = new CrMSIDataIO( basePath, fileNames, massChannels, numRows, dataType );
  CubeNumRows = numRows[0];
  cubes = new CrMSIDataIO::DataCube*[numOfThreadsDouble];
  iCube = new int[numOfThreadsDouble];
  bDataReady = new bool[numOfThreadsDouble];
  bRunningThread = new bool[numOfThreadsDouble];
  tworkers = new boost::thread[numOfThreadsDouble]; //There will be double of thread objects than the actually running threads
  bDataOverWrite = overWriteRamdisk;
  
}

ThreadingMsiProc::~ThreadingMsiProc()
{
  delete[] cubes;
  delete[] iCube;
  delete[] bDataReady;
  delete[] bRunningThread;
  delete[] tworkers; 
  delete ioObj;
}

void ThreadingMsiProc::runMSIProcessingCpp()
{
  //Initialize processing data cube
  life_end = false;
  for( int i = 0; i < numOfThreadsDouble; i++)
  {
    iCube[i] = -1; //-1 is not cube assigned to worker thread
    bDataReady[i] = false;
    bRunningThread[i] = false;
  }
  
  int nextCube = 0; //Point to the next datacube to load
  int runningThreads = 0; //Total number of running threads
  bool end_of_program = false;
  while( !end_of_program ) 
  {
    //Load data for future threads and start working threads
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      if(iCube[iThread] == -1 && nextCube < ioObj->getNumberOfCubes()) //No cube assigned then no thread running in this slot
      {
        Rcpp::Rcout<<"Processing cube "<<(nextCube + 1)<<" of "<<ioObj->getNumberOfCubes()<<"\n";
        iCube[iThread] = nextCube;
        nextCube++;
        cubes[iThread] = ioObj->loadDataCube(iCube[iThread]);
        if(2*runningThreads < numOfThreadsDouble)
        {
          bRunningThread[iThread] = true;
          tworkers[iThread] = boost::thread(boost::bind(&ThreadingMsiProc::ProcessingThread, this, iThread)); //Start Thread 
          runningThreads++;
        }
      }
    }
    
    //Wait for thread ends
    WaitForSomeThreadEnd();
    
    mtx.lock(); //Any thread locked to this mutex is actually waiting to save data
    
    //Update number of running threads
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      runningThreads = 0;
      if(bRunningThread[iThread])
      {
        runningThreads++;
      }
    }
    
    //Start thread that have data loaded ready to be processed
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      if(iCube[iThread] != -1 && !bRunningThread[iThread] && !bDataReady[iThread] && 2*runningThreads < numOfThreadsDouble)
      {
        bRunningThread[iThread] = true;
        tworkers[iThread] = boost::thread(boost::bind(&ThreadingMsiProc::ProcessingThread, this, iThread)); //Start Thread 
        runningThreads++;
      }      
    }
    
    //Save data to RSession and free thread slots
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      if( bDataReady[iThread] )
      {
        if(bDataOverWrite)
        {
          ioObj->storeDataCube(iCube[iThread], cubes[iThread]); //Overwrite datacube on HDD
        }
        ioObj->freeDataCube(cubes[iThread]); //Free datacube memory
        iCube[iThread] = -1; //Mark thread as stopped
        bDataReady[iThread] = false; //Reset data ready state;
        tworkers[iThread].join();
      }
    }
    
    //Check end condition
    if( nextCube >= ioObj->getNumberOfCubes() )
    {
      end_of_program = true;
      for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
      {
        end_of_program &= (iCube[iThread] == -1);
      }
    }
    mtx.unlock();
  }
  Rcpp::Rcout<<"Multi thread processing complete\n";
}

void ThreadingMsiProc::ProcessingThread( int threadSlot )
{
  //Call the processing function for this thread
  ProcessingFunction(threadSlot); 
  
  //Save data to R Session and Store the new state of total processed cubes
  mtx.lock();
  bDataReady[threadSlot] = true;
  bRunningThread[threadSlot] = false;
  mtx.unlock();
  
  //Notify Main thread
  boost::unique_lock<boost::mutex> lock(life_end_mutex);
  if(!life_end) //Notify only if it haven been done already by another thread
  {
    life_end = true;
    life_end_cond.notify_all();
  }
}

void ThreadingMsiProc::WaitForSomeThreadEnd()
{
  boost::unique_lock<boost::mutex> lock(life_end_mutex);
  while (!life_end)
  {
    life_end_cond.wait(lock);
  }
  life_end =false; //Reset it for the next thread
}

void ThreadingMsiProc::ProcessingFunction(int threadSlot)
{
  Rcpp::Rcout<<"Function ThreadingMsiProc::ProcessingFunction(int threadSlot) has been called from thread slot: "<<threadSlot<<"\n";
  Rcpp::Rcout<<"The base class ThreadingMsiProc can not be used directly and must be derived reimplementing the ProcessingFunction()\n";
}