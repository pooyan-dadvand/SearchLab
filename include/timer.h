/* -*- c++ -*- */
#pragma once

#include <ctime>
#include "omp.h"

double GetCurrentTime() {
#ifndef _OPENMP
  return static_cast< double >( std::clock()) / static_cast< double >( CLOCKS_PER_SEC );
#else
  return omp_get_wtime();
#endif
}
