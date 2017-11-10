#pragma once

#include <string>
#include <iostream>

#include "mpi.h"
#include "timer.h"

void printMaxTime(double a, double b, int mpi_rank, int mpi_size, std::string msg) {
  double local = b - a;
  double min, max;
  
  MPI_Allreduce(&local, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  double mid = (max + min) / 2.0;

  int local_up = ((b - a) >= mid);
  int local_dw = ((b - a) <= mid);

  int count_up, count_dw;

  MPI_Reduce(&local_up, &count_up, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_dw, &count_dw, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(!mpi_rank) {
    std::cout << "Min (" << mpi_rank << ") " << msg << ":\t" << min << " -- " << count_dw << std::endl;
    std::cout << "Max (" << mpi_rank << ") " << msg << ":\t" << max << " -- " << count_up << std::endl;
  }
}