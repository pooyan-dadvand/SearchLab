/* -*- c++ -*- */
#pragma once

#include <iostream>
#include <vector>

#include "point.h"
#include "timer.h"

#include "points_bins.h"
#include "points_bins_hash.h"

namespace PointsNew {

  template < class BinsType >
    void RunTests( char const Title[], std::vector< Point > &points_vector, const Point &search_point,
		   double radius, std::size_t numsearch, std::size_t numsearch_nearest,
		   const std::size_t GridSize[ 3], bool print_statistics ) {
    static bool first_time = true;
    double t0 = GetCurrentTime();
    BinsType bins( points_vector.begin(), points_vector.end(), GridSize );
    double t1 = GetCurrentTime();

    if ( print_statistics) {
      bins.PrintStatistics();
    }

    // Header of timmings
    if ( first_time ) {
      std::locale prev_loc = std::cout.getloc();
      std::cout.imbue( std::locale( "" ) ); // for thousand separators ...
      std::cout << "Doing " << numsearch << " searches in radius and " << std::endl;
      std::cout << "doing " << numsearch_nearest << " nearest searches." << std::endl;
      std::cout.imbue( prev_loc ); // restore previous locale, i.e. without thousand separators
      std::cout << "SS = Title\tGEN = generation\tSIR = search in radius\tSN = search nearest" << std::endl;
      std::cout << "OMP = OpenMP\tSER = serial\tNOFR = # points in SIR\tNP = nearest Point" << std::endl;
      std::cout << "SS\t\tGEN\tSIROMP\tSIRSER\tSNPOMP\tSNPSER\tNOFR\tNP" << std::endl;
      first_time = false;
    }
    // Title & GENeration
    std::cout << Title // "Points Bin"
	      << "\t" << t1 - t0 << "\t";

    std::vector< PointsBins< Point >::ResultType > results;
    PointsBins< Point >::ResultType nearest_point_result;

    // Search In Radius OpenMP
    t0 = GetCurrentTime();
#pragma omp parallel for firstprivate( results )
    for ( std::size_t i = 0; i < numsearch; i++ ) {
      results.clear();
      bins.SearchInRadius( search_point, radius, results );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    // Search In Radius SERial
    t0 = GetCurrentTime();
    for ( std::size_t i = 0; i < numsearch; i++ ) {
      results.clear();
      bins.SearchInRadius( search_point, radius, results );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    // Search Nearest Parallel OpenMP
    t0 = GetCurrentTime();
#pragma omp parallel for firstprivate( nearest_point_result )
    for ( std::size_t i = 0; i < numsearch_nearest; i++ ) {
      nearest_point_result = bins.SearchNearest( search_point );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    // Search Nearest Parallel? SERial + NOFR? and Neares Point 
    t0 = GetCurrentTime();
    for ( std::size_t i = 0; i < numsearch_nearest; i++ ) {
      nearest_point_result = bins.SearchNearest( search_point );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t" << results.size() << "\t" << *( nearest_point_result.Get());
    std::cout << std::endl;
  }
}
