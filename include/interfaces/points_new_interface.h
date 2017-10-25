/* -*- c++ -*- */
#pragma once

#include <iostream>
#include <vector>

#include "point.h"
#include "timer.h"

#include "points_bins.h"

namespace PointsNew {

  template < class BinsType >
    void RunTests( char const Title[], std::vector< Point > &points_vector, const Point &search_point,
		   double radius, std::size_t numsearch, std::size_t numsearch_nearest ) {
    static bool first_time = true;
    double t0 = GetCurrentTime();
    BinsType bins( points_vector.begin(), points_vector.end() );
    double t1 = GetCurrentTime();

    bins.PrintStatistics();

    // Header of timmings
    if ( first_time ) {
      std::cout << "SS\t\tGEN\tSIROMP\tSIRSER\tSNPOMP\tSNPSER\tNOFR\tNP" << std::endl;
      first_time = false;
    }
    std::cout << Title // "Points Bin"
	      << "\t" << t1 - t0 << "\t";

    std::vector< PointsBins< Point >::ResultType > results;
    PointsBins< Point >::ResultType nearest_point_result;

    t0 = GetCurrentTime();
#pragma omp parallel for firstprivate( results )
    for ( std::size_t i = 0; i < numsearch; i++ ) {
      results.clear();
      bins.SearchInRadius( search_point, radius, results );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    t0 = GetCurrentTime();
    for ( std::size_t i = 0; i < numsearch; i++ ) {
      results.clear();
      bins.SearchInRadius( search_point, radius, results );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    t0 = GetCurrentTime();

#pragma omp parallel for firstprivate( nearest_point_result )
    for ( std::size_t i = 0; i < numsearch_nearest; i++ ) {
      nearest_point_result = bins.SearchNearest( search_point );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    t0 = GetCurrentTime();
    for ( std::size_t i = 0; i < numsearch_nearest; i++ ) {
      nearest_point_result = bins.SearchNearest( search_point );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t" << results.size() << "\t" << *nearest_point_result.Get();
    std::cout << std::endl;
  }
}
