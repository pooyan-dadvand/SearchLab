#pragma once

#include <ctime>

#include "omp.h"

double GetCurrentTime() {
#ifndef _OPENMP
	return std::clock() / static_cast<double>(CLOCKS_PER_SEC);
#else
	return omp_get_wtime();
#endif
}
