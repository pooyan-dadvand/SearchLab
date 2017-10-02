#pragma once

// System includes>
#include <string>
#include <iostream>
#include <algorithm>

#include <assert.h> 

// OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Mpi
#include "mpi.h"
#include "partition_bins.h"
#include "global_pointer.h"

// Timer
#include "timer_mpi.h"

/**
 * Mpi version of the Bins.
 */
template<class TObjectBins, class TPartitionBins = PartitionBins<typename TObjectBins::ObjectType>>
class ParallelBins {
  static constexpr int Dimension = 3;

public:

  // Use the types from the hosted bins
  using InternalPointType = typename TObjectBins::InternalPointType;
	using PointerType = typename TObjectBins::PointerType;
	using ResultType = typename TObjectBins::ResultType;
  using ObjectType = typename TObjectBins::ObjectType;
  using ResultArray = std::vector<ResultType>;

  /** Constructor using the objects.
   * Creates a bins by calculating the cell size based on the objects
   * from 'ObjectsBegin' to 'ObjectsEnd' and fills it with such objects.
   * @param ObjectsBegin Iterator pointing to the first object to add to the bins
   * @param ObjectsEnd   Iterator pointing to the last object to add to the bins
   */
  template<typename TIteratorType>
  ParallelBins(TIteratorType const& LocalObjectsBegin, TIteratorType const& LocalObjectsEnd) {
    InitializeMpi();

    GenerateLocalBins(LocalObjectsBegin, LocalObjectsEnd);
    GeneratePartitionBins(LocalObjectsBegin, LocalObjectsEnd, mpi_size * 10);
  }

  /** Assignment operator.
   * Assignment operator.
   * @param rOther reference object
   */
  ParallelBins<TObjectBins, TPartitionBins> & operator=(const ParallelBins<TObjectBins, TPartitionBins> & rOther) {
    mPartitionBins = rOther.mPartitionBins;
    mObjectBins = rOther.mObjectBins;

    return *this;
  }

  /** Copy constructor.
   * Copy constructor.
   * @param rOther reference object
   */
  ParallelBins(const ParallelBins& rOther) {
    *this = rOther;
  }

  /** Default destructor
   * Default destructor
   */
  virtual ~ParallelBins() {
    delete mPartitionBins;
    delete mObjectBins;
  }

  /** Search in radius
   * Search in radius
   * @param  ObjectsBegin Iterator to the first object to search
   * @param  ObjectsEnd   Iterator to the last object to search
   * @param  Radius       Radius of the search
   * @param  ResultsArray Array of ResultsArray ( 1 resultsArray per object searched )
   */
  template<typename TIteratorType>
  void SearchInRadius(
      const TIteratorType & ObjectsBegin,
      const TIteratorType & ObjectsEnd,
      double Radius,
      std::vector<ResultArray> & ResultsArray) {

    std::size_t NumberOfObjects = ObjectsEnd - ObjectsBegin;

    double t0, t1, ti;

    // Not sure what happens here if a partition ends up having 0 points/objects/elements
    assert(NumberOfObjects > 0);

    std::vector<std::vector<ObjectType>>  SendObjects(mpi_size, std::vector<ObjectType>(0));
    std::vector<std::vector<ObjectType>>  RecvObjects(mpi_size, std::vector<ObjectType>(0));
    std::vector<std::vector<ResultArray>> SendResults(mpi_size, std::vector<ResultArray>(0));
    std::vector<std::vector<ResultArray>> RecvResults(mpi_size, std::vector<ResultArray>(0));

    std::vector<std::vector<double>> SendRadius(mpi_size, std::vector<double>(0));
    std::vector<std::vector<double>> RecvRadius(mpi_size, std::vector<double>(0));

    std::vector<std::vector<int>> SentObjectsIds(mpi_size, std::vector<int>(0));

    MPI_Request * sendReq = new MPI_Request[mpi_size];
    MPI_Request * recvReq = new MPI_Request[mpi_size];

    MPI_Status * sendStat = new MPI_Status[mpi_size];
    MPI_Status * recvStat = new MPI_Status[mpi_size];

    // Calculate the remote partitions where we need to search for each point and execute the transfer in background
    mPartitionBins->cleanDebugStorage();
    t0 = GetCurrentTime();
    SearchPartitions(ObjectsBegin, NumberOfObjects, Radius, SendObjects, SendRadius, SentObjectsIds);
    t1 = GetCurrentTime();
    mPartitionBins->printDebugStorage();

    printMaxTime(t0, t1, mpi_rank, mpi_size, "Partition search time");

    t0 = GetCurrentTime();
    SendTransferPoints(SendObjects, sendReq, recvReq);
    t1 = GetCurrentTime();

    printMaxTime(t0, t1, mpi_rank, mpi_size, "Send spawn time");

    t0 = GetCurrentTime();
    SearchInRadiusLocal(ObjectsBegin, NumberOfObjects, Radius, ResultsArray);
    t1 = GetCurrentTime();

    printMaxTime(t0, t1, mpi_rank, mpi_size, "Loacal search time");

    t0 = GetCurrentTime();
    // Recv the points after the local search has finished and execute the remote search (all points should be here by now)
    for(int p = 0; p < mpi_size; p++) {
      RecvTransferPoints(RecvObjects, sendReq, recvReq, sendStat, recvStat, p);
    }
    t1 = GetCurrentTime();

    printMaxTime(t0, t1, mpi_rank, mpi_size, "Recv gather time");

    t0 = GetCurrentTime();
    for(int p = 0; p < mpi_size; p++) {
      if(p != mpi_rank && RecvObjects[p].size() != 0) {
        auto ObjectsRemoteBeg = RecvObjects[p].begin();
        auto ObjectsRemoteEnd = RecvObjects[p].end();
        std::size_t NumberOfRecvObjects = ObjectsRemoteEnd - ObjectsRemoteBeg;

        SendResults[p].resize(RecvObjects[p].size());

        SearchInRadiusLocal(ObjectsRemoteBeg, NumberOfRecvObjects, Radius, SendResults[p]);
      }
    }
    t1 = GetCurrentTime();

    printMaxTime(t0, t1, mpi_rank, mpi_size, "Remote search time");
    
    // Send the results (global pointers) back to the propper process and merge with the localss
    t0 = GetCurrentTime();
    TransferResults(SendObjects, SendResults, RecvResults);
    t1 = GetCurrentTime();

    printMaxTime(t0, t1, mpi_rank, mpi_size, "Transfer res time");

    t0 = GetCurrentTime();
    AssembleResults(SentObjectsIds, ResultsArray, RecvResults);
    t1 = GetCurrentTime();

    printMaxTime(t0, t1, mpi_rank, mpi_size, "Aseemble results");

    std::size_t accumNumResults = 0;
    for(std::size_t i = 0; i < NumberOfObjects; i++) {
      accumNumResults += ResultsArray[i].size();
      // std::cout << "(" << mpi_rank << ") Point (" << (*(ObjectsBegin + i))[0] << "," << (*(ObjectsBegin + i))[1] << "," << (*(ObjectsBegin + i))[2] << ") results: " << ResultsArray[i].size() << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(!mpi_rank) {
      std::cout << "Finish: " << ResultsArray[NumberOfObjects/2].size() << " " << accumNumResults << std::endl;
    }

    delete[] sendReq;
    delete[] recvReq;

    delete[] sendStat;
    delete[] recvStat;
  }

private:

  /**
   * Default constructor its private and should never be called
   */
  ParallelBins() {
    std::cerr << "This should never happen" << std::endl;
    abort();
  }

  /**
   * Initializes mpi_rank and mpi_size for the given mpi communicator.
   * If no mpi communicator is provided MPI_COMM_WORLD is used as default.
   */
  void InitializeMpi(MPI_Comm mpiComm = MPI_COMM_WORLD) {
    MPI_Comm_rank(mpiComm, &mpi_rank);
    MPI_Comm_size(mpiComm, &mpi_size);
  }

  /**
   * Generate the bins for the local objects
   * @param  ObjectsBegin Iterator to the first local object
   * @param  ObjectsEnd   Iterator to the last local object
   */
  template<typename TIteratorType>
  void GenerateLocalBins(TIteratorType const& ObjectsBegin, TIteratorType const& ObjectsEnd) {
    mObjectBins = new TObjectBins(ObjectsBegin, ObjectsEnd);
  }

  /**
   * Generate the bins for the local objects
   * @param ObjectsBegin Iterator to the first local object
   * @param ObjectsEnd   Iterator to the last local object
   * @param PartitionBinsSize Orientative size of the bins
   */
  template<typename TIteratorType>
  void GeneratePartitionBins(TIteratorType const& ObjectsBegin, TIteratorType const& ObjectsEnd, std::size_t PartitionBinsSize) {
    mPartitionBins = new TPartitionBins(ObjectsBegin, ObjectsEnd, PartitionBinsSize);
  }

  /** SearchPartitions
   * SearchInRadiusLocal
   * @param  ObjectsBegin       [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  template<typename TIteratorType>
  void SearchPartitions(
      const TIteratorType & ObjectsBegin,
      const std::size_t & NumberOfObjects,
      double const & Radius,
      std::vector<std::vector<ObjectType>> & SendObjects,
      std::vector<std::vector<double>> & SendRadius,
      std::vector<std::vector<int>> & SendObjectsIds) {

    // int TotalToSend = 0;
    // int GlobalToSend = 0;

    // int TotalNumObjects = NumberOfObjects;
    // int GlobalNumObjects = 0;

    double t0, t1, tx, tz;
    double accum_a = 0.0;
    double accum_b = 0.0;

    tx = GetCurrentTime();
    for(std::size_t i = 0; i < NumberOfObjects; i++) {
      auto ObjectItr = ObjectsBegin + i;

      t0 = GetCurrentTime();
      std::vector<typename TPartitionBins::ResultType> partitionList;
      mPartitionBins->SearchInRadius(*ObjectItr, Radius, partitionList);
      t1 = GetCurrentTime();

      accum_a += (t1 - t0);

      t0 = GetCurrentTime();
      for(std::size_t j = 0; j < partitionList.size(); j++) {
        int part = partitionList[j];
        if(part != mpi_rank && part != -1){
          SendObjects[part].push_back(*ObjectItr);
          SendRadius[part].push_back(Radius);
          SendObjectsIds[part].push_back(i);
          // TotalToSend++;
        }
      }
      t1 = GetCurrentTime();
      accum_b += (t1 - t0);
    }
    tz = GetCurrentTime();

    printMaxTime(0.0, accum_a, mpi_rank, mpi_size, "Search? ");
    printMaxTime(0.0, accum_b, mpi_rank, mpi_size, "Arrays? ");
    printMaxTime(tx, tz, mpi_rank, mpi_size, "Total? ");

    // MPI_Allreduce(&TotalToSend, &GlobalToSend, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // MPI_Allreduce(&TotalNumObjects, &GlobalNumObjects, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // std::cout << "(" << mpi_rank << ") " << "Sending " << TotalToSend << " unic objects out of " << TotalNumObjects << std::endl;
    // std::cout << "(" << mpi_rank << ") " << "Sending " << GlobalToSend << " unic objects out of " << GlobalNumObjects << std::endl;
  }

  /** SearchInRadiusLocal
   * SearchInRadiusLocal
   * @param  ThisObject         [description]
   * @param  Radius             [description]
   * @param  Results            [description]
   * @param  MaxNumberOfResults [description]
   * @return                    [description]
   */
  template<typename TIteratorType>
  void SearchInRadiusLocal(
      const TIteratorType & ObjectsBegin,
      const std::size_t & NumberOfObjects,
      double const & Radius,
      std::vector<ResultArray> & Results) {

    auto ResultsBegin = Results.begin();

    for(std::size_t i = 0; i < NumberOfObjects; i++) {
      auto ObjectItr  = ObjectsBegin + i;
      auto ResultsItr = ResultsBegin + i;
      mObjectBins->SearchInRadius(*ObjectItr, Radius, *ResultsItr);
      // std::cout << "Point:" << " (" << (*ObjectItr)[0] << "," << (*ObjectItr)[1] << "," << (*ObjectItr)[2] << ") " << "Found: " <<  (*ResultsItr).size() << " Results." << std::endl;
    }
  }

  void SendTransferPoints(const std::vector<std::vector<ObjectType>> & sendObjects,
      MPI_Request * sendReq, MPI_Request * recvReq) {

    sendSize = new int[mpi_size];
    recvSize = new int[mpi_size];

    sendBuffers = std::vector<double *>(mpi_size, nullptr);
    recvBuffers = std::vector<double *>(mpi_size, nullptr);

    
    for(int p = 0; p < mpi_size; p++) {
      sendSize[p] = sendObjects[p].size();
    }

    MPI_Alltoall(sendSize, 1, MPI_INT, recvSize, 1, MPI_INT, MPI_COMM_WORLD);
    // std::cout << "(" << mpi_rank << ") sendrecvSize: ";
    // for(int p = 0; p < mpi_size; p++) {
    //   std::cout << "[" << p << "] " << sendSize[p] << "|" << sendSize[p];
    // }
    // std::cout << std::endl;

    for(int p = 0; p < mpi_size; p++) {
      sendBuffers[p] = new double[sendSize[p] * 3];
      recvBuffers[p] = new double[recvSize[p] * 3];

      for(std::size_t i = 0; i < sendObjects[p].size(); i++) {
        sendBuffers[p][i * 3 + 0] = sendObjects[p][i][0];
        sendBuffers[p][i * 3 + 1] = sendObjects[p][i][1];
        sendBuffers[p][i * 3 + 2] = sendObjects[p][i][2];
      }

      MPI_Isend(sendBuffers[p], 3 * sendSize[p], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &sendReq[p]);
      MPI_Irecv(recvBuffers[p], 3 * recvSize[p], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &recvReq[p]);
    }
  }

  void RecvTransferPoints(std::vector<std::vector<ObjectType>> & recvObjects,
      MPI_Request * sendReq, MPI_Request * recvReq, MPI_Status * sendStat, MPI_Status * recvStat, int p) {
    // Wait until everything has been sent and recv
    MPI_Wait(&sendReq[p], &sendStat[p]);
    MPI_Wait(&recvReq[p], &recvStat[p]);

    recvObjects[p].resize(recvSize[p]);

    for(std::size_t j = 0; j < recvObjects[p].size(); j++) {
      recvObjects[p][j][0] = recvBuffers[p][j * 3 + 0];
      recvObjects[p][j][1] = recvBuffers[p][j * 3 + 1];
      recvObjects[p][j][2] = recvBuffers[p][j * 3 + 2];
      // std::cout << "(" << mpi_rank << ") " << "got point " << "(" << (*(recvPoints[i][j]))[0] << "," << (*(recvPoints[i][j]))[1] << "," << (*(recvPoints[i][j]))[2] << ") from proc: " << i << std::endl; 
    }

    delete[] sendBuffers[p];
    delete[] recvBuffers[p];
  }

  /** Sends the results back to the right process.
   * Sends the results back to the right process.
   * @param  sendObjects Objects that were sent to this processor.
   * @param  sendResults List of results to send to anoher processor.
   * @param  recvResults List of results that have been sent to this processor.
   */
  void TransferResults(const std::vector<std::vector<ObjectType>> & sendObjects, std::vector<std::vector<ResultArray>> & sendResults, std::vector<std::vector<ResultArray>> & recvResults) {
    
    MPI_Request * sendReq = new MPI_Request[mpi_size];
    MPI_Request * recvReq = new MPI_Request[mpi_size];

    MPI_Status * sendStat = new MPI_Status[mpi_size];
    MPI_Status * recvStat = new MPI_Status[mpi_size];

    // Get the number of the results for every point for every processor and communicate them
    std::vector<int *> sendNumResults(mpi_size, nullptr);
    std::vector<int *> recvNumResults(mpi_size, nullptr);

    // Resize the results to send
    for(int p = 0; p < mpi_size; p++) {
      sendNumResults[p] = new int[sendResults[p].size()];
      for(std::size_t i = 0; i < sendResults[p].size(); i++) {
        sendNumResults[p][i] = sendResults[p][i].size();
      }
    }

    // Resize the results to recieve ( we will recieve from N as much results as points we sent to the proc N)
    for(int p = 0; p < mpi_size; p++) {
      recvNumResults[p] = new int[sendObjects[p].size()];
    }

    // Spawn the messages
    for(int p = 0; p < mpi_size; p++) {
      MPI_Isend(sendNumResults[p], sendResults[p].size(), MPI_INT, p, 0, MPI_COMM_WORLD, &sendReq[p]);
      MPI_Irecv(recvNumResults[p], sendObjects[p].size(), MPI_INT, p, 0, MPI_COMM_WORLD, &recvReq[p]);
    }

    // Get the size of a global pointer (in char)
    std::size_t gp_size = sizeof(GlobalPointer<ObjectType>);

    std::vector<char *> sendResultData(mpi_size, nullptr);
    std::vector<char *> recvResultData(mpi_size, nullptr);

    for(int p = 0; p < mpi_size; p++) {
      MPI_Wait(&sendReq[p], &sendStat[p]);
      MPI_Wait(&recvReq[p], &recvStat[p]);
    }

    // At this point we know the number of results we will recieve for each point from every processor.
    // Prepare the buffers to send and recieve that data.

    // Results are sent PER process and not per point, this is done to reduce the number of messages from N * P to P
    // While the size will remain the same the number of messages is much lower, so the network does not saturate
    for(int p = 0; p < mpi_size; p++) {
      std::size_t accumSendRes = 0;
      std::size_t accumRecvRes = 0;

      for(std::size_t i = 0; i < sendResults[p].size(); i++) {
        accumSendRes += sendNumResults[p][i];
      }

      for(std::size_t i = 0; i < sendObjects[p].size(); i++) {
        accumRecvRes += recvNumResults[p][i];
      }

      // Only send if there is something to send ( avoid having P*P messages )
      if(accumSendRes) {
        sendResultData[p] = new char[gp_size * accumSendRes];
        for(std::size_t i = 0; i < sendResults[p].size(); i++) {        // The particles
          for(std::size_t j = 0; j < sendResults[p][i].size(); j++) {   // The results
            sendResults[p][i][j].Get().Save(&sendResultData[p][j * gp_size]);
          }
        }
        MPI_Isend(sendResultData[p], gp_size * accumSendRes, MPI_CHAR, p, 0, MPI_COMM_WORLD, &sendReq[p]);
      }

      // Only recv if there is something to recv ( avoid having P*P messages )
      if(accumRecvRes) {
        recvResultData[p] = new char[gp_size * accumRecvRes];
        MPI_Irecv(recvResultData[p], gp_size * accumRecvRes, MPI_CHAR, p, 0, MPI_COMM_WORLD, &recvReq[p]);
      }
    }

    // Wait until all communication has been done
    for(int p = 0; p < mpi_size; p++) {
      MPI_Wait(&sendReq[p], &sendStat[p]);
      MPI_Wait(&recvReq[p], &recvStat[p]);
    }

    for(int p = 0; p < mpi_size; p++) {
      recvResults[p].resize(sendObjects[p].size());
      for(std::size_t i = 0; i < recvResults[p].size(); i++) {        // The particles
        recvResults[p][i].resize(recvNumResults[p][i]);
        for(std::size_t j = 0; j < recvResults[p][i].size(); j++) {   // The results
          recvResults[p][i][j].Get().Load(&recvResultData[p][j * gp_size]);
        }
      }
    }

    // Free the memory
    for(int p = 0; p < mpi_size; p++) {
      delete[] sendNumResults[p];
      delete[] recvNumResults[p];

      if(sendResultData[p] != nullptr) {
        delete[] sendResultData[p];
      }
      if(recvResultData[p] != nullptr) {
        delete[] recvResultData[p];
      }
    }

    delete[] sendReq;
    delete[] recvReq;
    delete[] sendStat;
    delete[] recvStat;
  }

  /** Merges remote and local results.
   * Merges remote and local results. Results will be added to 'Results'
   * @param  SentObjectsMap Map of objects sent to processes.
   * @param  Results        List of local results for this process search points.
   * @param  recvResults    List of remote results for this process search points.
   */
  void AssembleResults(std::vector<std::vector<int>> & SendObjectsIds, std::vector<ResultArray> & Results, std::vector<std::vector<ResultArray>> & recvResults) {
    for(int p = 0; p < mpi_size; p++) {
      auto sentToProc = SendObjectsIds[p].size();
      if(sentToProc) {
        std::size_t recvResIndex = 0;
        for(std::size_t i = 0; i < sentToProc; i++) { // Iterate over the results arrays
          auto numRemResults = recvResults[p][recvResIndex].size();
          if(numRemResults > 0) { // It still can be 0!!!
            Results[SendObjectsIds[p][i]].reserve(Results[SendObjectsIds[p][i]].size()+numRemResults);
            Results[SendObjectsIds[p][i]].insert(Results[SendObjectsIds[p][i]].end(),recvResults[p][recvResIndex].begin(), recvResults[p][recvResIndex].end());
          }
          recvResIndex++;
        }
      }
    }
  }

  /** Calculates the squared distance.
   * Calculates the squared distance between two points.
   */
  double Distance2(ObjectType const& FirstPoint, ObjectType const& SecondPoint) {
		double result = double();
		for (int i = 0; i < Dimension; i++) {
			auto distance_i = FirstPoint[i] - SecondPoint[i];
			result += distance_i * distance_i;
		}
		return result;
	}
 
  /** Mpi related variables
   * Mpi related variables
   * @mpi_rank: id of the current process for the given mpi communicator
   * @mpi_size: number of processes for the given mpi communicator
   */
  int mpi_rank;
  int mpi_size;

  /**
   * Bins used in the BinsMpi.
   * @mPartitionBin: Stores where the partitions are located in the space
   * @mObjectBin: Stores where the objects of the partition/s associated to this process are located in the space
   */
  TPartitionBins  * mPartitionBins;
  TObjectBins     * mObjectBins;
  std::size_t       mObjectsSize;
  
  /**
   * Async transfer stuff
   **/
  int * sendSize;
  int * recvSize;

  std::vector<double *> sendBuffers;
  std::vector<double *> recvBuffers;
};

/// input stream function
template<class TObjectConfigure, class TPartConfigure>
inline std::istream& operator >> (
    std::istream& rIStream,
    ParallelBins<TObjectConfigure>& rThis) {

  return rIStream;
}


/// output stream function
template<class TObjectConfigure, class TPartConfigure>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ParallelBins<TObjectConfigure, TPartConfigure> & rThis) {

  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
