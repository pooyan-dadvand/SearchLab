#include <iostream>
#include <vector>

#include "timer.h"
#include "point.h"
#include "mpi.h"

#include "points_bins.h"
#include "global_pointer.h"

namespace PointsNew {

template<class BinsType>
void RunTests(char const Title[], std::vector<Point> & points_vector, const Point & search_point, double radius, std::size_t numsearch, std::size_t numsearch_nearest) {

  int mpi_rank, mpi_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  std::size_t numGlobPoints = points_vector.end() - points_vector.begin();

  int numBigProc = numGlobPoints % mpi_size;                                // Procs with haigh balance
  int numLowProc = mpi_size - numBigProc;                                   // Procs with low balance
  int sizeBigProc = std::ceil((double)numGlobPoints / (double)mpi_size);
  int sizeLowProc =  std::floor((double)numGlobPoints / (double)mpi_size);

  std::size_t numProcPoints = mpi_rank < numBigProc ? sizeBigProc : sizeLowProc;
  std::size_t pointOffset = (mpi_rank < numBigProc) * (mpi_rank * sizeBigProc) + (mpi_rank >= numBigProc) * ((numBigProc * sizeBigProc) + ((mpi_rank - numBigProc) * sizeLowProc));

  auto pointsPartBeg = points_vector.begin() + pointOffset;
  auto pointsPartEnd = pointsPartBeg + numProcPoints;

  // std::cout << "numBigProc: " << numBigProc << std::endl;
  // std::cout << "numLowProc: " << numLowProc << std::endl;
  // std::cout << "sizeBigProc: " << sizeBigProc << std::endl;
  // std::cout << "sizeLowProc: " << sizeLowProc << std::endl;
  // std::cout << "numProcPoints: " << numProcPoints << std::endl;
  // std::cout << "pointOffset: " << pointOffset << std::endl;
  // std::cout << "=============" << std::endl;

  // Synchronize all the processes before timing
  MPI_Barrier(MPI_COMM_WORLD);

  double t0 = GetCurrentTime();
  BinsType bins(pointsPartBeg, pointsPartEnd);
  double t1 = GetCurrentTime();

  std::cout << t1 - t0 << "\t";

  std::vector<std::vector<PointsBins<Point>::ResultType>> resultsArray(numProcPoints, std::vector<PointsBins<Point>::ResultType>(0));

  PointsBins<Point>::ResultType nearest_point_result;

  t0 = GetCurrentTime();
  bins.SearchInRadius(pointsPartBeg, pointsPartEnd, radius, resultsArray);
  t1 = GetCurrentTime();

  std::cout << t1 - t0 << "\t";

  // t0 = GetCurrentTime();
  // for (std::size_t i = 0; i < numsearch; i++) {
  //   results.clear();
  //   bins.SearchInRadius(search_point, radius, results);
  // }
  // t1 = GetCurrentTime();
  // std::cout << t1 - t0 << "\t";

  // t0 = GetCurrentTime();

  // #pragma omp parallel for firstprivate(nearest_point_result)
  // for (std::size_t i = 0; i < numsearch_nearest; i++) {
  //   nearest_point_result = bins.SearchNearest(search_point);
  // }
  // t1 = GetCurrentTime();
  // std::cout << t1 - t0 << "\t";

  // t0 = GetCurrentTime();
  // for (std::size_t i = 0; i < numsearch_nearest; i++) {
  //   nearest_point_result = bins.SearchNearest(search_point);
  // }
  // t1 = GetCurrentTime();
  // std::cout << t1 - t0 << "\t" << results.size() << "\t" << *nearest_point_result.Get();
  std::cout << std::endl;
}

}
