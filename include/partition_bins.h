#pragma once

#include <array>  
#include <vector>
#include <limits>
#include <cmath>

#include "mpi.h"

#include "partition_bins_cells_container.h"
#include "spatial_search_result.h"

#include "timer_mpi.h"

template <typename TObjectType>
class PartitionBins {
  static constexpr int Dimension = 3;
public:
  using InternalPointType = std::array<double, Dimension>;
  using ResultType = int;
  
  int mpi_size;
  int mpi_rank;

  template<typename TIteratorType>
  PartitionBins(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd, std::size_t numGlobPoints) : mCells(PointsBegin, PointsEnd, numGlobPoints) {
    mNumberOfPoints = std::distance(PointsBegin, PointsEnd);
    
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (mNumberOfPoints == 0) {
      mpPartitions = nullptr;
      return;
    }
    
    // At this point mCells contains an array of indices which by construction are 0 or 1. If we reduce these
    // arrays we "merge" the virtual vector poiting at the partitions and we can calculate its size
    std::size_t numberOfCells = mCells.GetTotalNumberOfCells();

    int * sendOffsets   = new int[numberOfCells];
    int * recvOffsets   = new int[numberOfCells * mpi_size];
    int * globalOffsets = new int[numberOfCells];

    std::size_t sendReqIdx = 0;
    std::size_t recvReqIdx = 0;

    MPI_Request * sendReq = new MPI_Request[mpi_size];
    MPI_Request * redvReq = new MPI_Request[mpi_size];

    MPI_Status * sendStat = new MPI_Status[mpi_size];
    MPI_Status * recvStat = new MPI_Status[mpi_size];

    for(std::size_t i = 0; i < numberOfCells; i++) {
      sendOffsets[i] = mCells.GetCellBeginIndex(i);
      for(int j = 0; j < mpi_size; j++) {
        recvOffsets[j * numberOfCells + i] = 0;
      }
    }

    // Spawn the collectives
    for(int i = 0; i < mpi_size; i++) {
      int srcIdx = 0;
      int dstIdx = numberOfCells * i;
      
      MPI_Isend(&sendOffsets[srcIdx], numberOfCells, MPI_INT, i, 0, MPI_COMM_WORLD, &sendReq[sendReqIdx++]);
      MPI_Irecv(&recvOffsets[dstIdx], numberOfCells, MPI_INT, i, 0, MPI_COMM_WORLD, &redvReq[recvReqIdx++]);
    }

    // Wait for all send
    for(int p = 0; p < mpi_size; p++) {
      MPI_Wait(&sendReq[p], &sendStat[p]);
      MPI_Wait(&redvReq[p], &recvStat[p]);
    }

    // Initialize the offset array
    for(std::size_t i = 0; i < numberOfCells; i++) {
      globalOffsets[i] = 0;
    }

    for(int p = 0; p < mpi_size; p++) {
      for(std::size_t i = 0; i < numberOfCells; i++) {
        globalOffsets[i] += recvOffsets[numberOfCells * p + i];
      }
    }

    // Update the bins
    for(std::size_t i = 0; i < numberOfCells; i++) {
      mCells.SetCellBeginIndex(i, globalOffsets[i]);
    }

    // Create the array with the partitions
    std::size_t NumPartsPerCell = globalOffsets[numberOfCells-1];
    
    mpPartitions = new int[NumPartsPerCell];

    int * sendPartitions = new int[NumPartsPerCell];
    int * recvPartitions = new int[NumPartsPerCell];

    for (std::size_t i = 0; i < NumPartsPerCell; i++) {
      mpPartitions[i] = -1;
    }

    AssignPartitionsToCells(PointsBegin, PointsEnd, recvOffsets, numberOfCells);

    for (std::size_t i = 0; i < NumPartsPerCell; i++) {
      sendPartitions[i] = mpPartitions[i];
    }

    MPI_Allreduce(sendPartitions, recvPartitions, NumPartsPerCell, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    for (std::size_t i = 0; i < NumPartsPerCell; i++) {
      mpPartitions[i] = recvPartitions[i];
    }

    delete[] sendOffsets;
    delete[] recvOffsets;
    delete[] globalOffsets;

    delete[] sendReq;
    delete[] redvReq;

    delete[] sendStat;
    delete[] recvStat;

    delete[] sendPartitions;
    delete[] recvPartitions;
  }

  void SearchInRadius(TObjectType const& ThePoint, double Radius, std::vector<ResultType>& rResults) {
    InternalPointType min_point;
    std::array<std::size_t, Dimension> length;

    for (int i = 0; i < Dimension; i++) {
      min_point[i] = ThePoint[i] - Radius;
      length[i] = mCells.CalculatePosition(ThePoint[i] + Radius, i) - mCells.CalculatePosition(ThePoint[i] - Radius, i) + 1;
    }
    auto min_cell = mCells.CalculateCellIndex(min_point);

    std::vector<int> usedMask(mpi_size, 0);

    for (std::size_t i_z = 0; i_z < length[2]; i_z++) {
      auto y_position = min_cell + i_z * mCells.GetNumberOfCells(0) * mCells.GetNumberOfCells(1);
      for (std::size_t i_y = 0; i_y < length[1]; i_y++) {
        std::size_t offset = mCells.GetCellBeginIndex(y_position);
        const std::size_t end_offset = mCells.GetCellBeginIndex(y_position + length[0]);
        int * p_part = mpPartitions + mCells.GetCellBeginIndex(y_position);
        for (; offset <end_offset; offset++) {
          if(!usedMask[*p_part] && *p_part != -1)
          { 
            usedMask[*p_part] = 1;
            rResults.push_back(*p_part);
          }
          p_part++;
        }
        y_position += mCells.GetNumberOfCells(0);
      }
    }
  }

private:
  std::size_t mNumberOfPoints;
  PartitionBinsCellsContainer mCells;
  int * mpPartitions;

  template<typename TIteratorType>
  void AssignPartitionsToCells(TIteratorType const& PointsBegin, TIteratorType const& PointsEnd, int * recvOffsets, std::size_t numberOfCells) {
    for (auto i_point = PointsBegin; i_point != PointsEnd; i_point++) {
      auto index = mCells.CalculateCellIndex(*i_point);
      std::size_t offset = mCells.GetCellBeginIndex(index);

      for(int p = 0; p < mpi_size; p++) {
        if((recvOffsets[numberOfCells * p + index + 1] - recvOffsets[numberOfCells * p + index]) == 1) {
          mpPartitions[offset++] = p;
        }
      }
    }
  }

  double Distance2(TObjectType const& FirstPoint, TObjectType const& SecondPoint) {
    double result = double();
    for (int i = 0; i < Dimension; i++) {
      auto distance_i = FirstPoint[i] - SecondPoint[i];
      result += distance_i * distance_i;
    }
    return result;
  }
};