/* -*- c++ -*- */
#pragma once

#include <iostream>
#include <vector>

#include "point.h"
#include "timer.h"

#include "spatial_containers/spatial_containers.h"

namespace PointsOld {

  template < class SearchStructureType, class PointType, class IteratorType, class DistanceIterator >
    void RunTests( char const Title[], IteratorType PBegin, IteratorType PEnd, IteratorType results0,
		   DistanceIterator distances0, std::size_t MaxResults, PointType *allPoints,
		   double radius0, std::size_t numsearch, std::size_t bucket_size ) {
    static bool first_time = true;
    double t0, t1;
    std::size_t n;

    std::size_t npoints = PEnd - PBegin;
    std::vector< std::size_t > numresArray( numsearch );
    PointType &point0 = allPoints[ 0 ];

    std::size_t numsearch_nearest = 100 * numsearch;

    t0 = GetCurrentTime();
    SearchStructureType nodes_tree( PBegin, PEnd );
    t1 = GetCurrentTime();

    // Header of timmings
    if ( first_time ) {
      std::cout << "SS\t\tGEN\tSIROMP\tSIRSER\tSNPOMP\tSNPSER\tNOFR\tNP" << std::endl;
      first_time = false;
    }
    std::cout << Title << "\t" << t1 - t0 << "\t";

    t0 = GetCurrentTime();
#pragma omp parallel for
    for ( std::size_t i = 0; i < numsearch; i++ ) {
      IteratorType resultsItr = results0;
      DistanceIterator distanceItr = distances0;
      n = nodes_tree.SearchInRadius( point0, radius0, resultsItr, distanceItr, MaxResults );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    t0 = GetCurrentTime();
    for ( std::size_t i = 0; i < numsearch; i++ ) {
      IteratorType resultsItr = results0;
      DistanceIterator distanceItr = distances0;
      n = nodes_tree.SearchInRadius( point0, radius0, resultsItr, distanceItr, MaxResults );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    PointType **PNearestArray = new PointType *[ numsearch ];
    std::vector< double > distancesArrayNearest( numsearch );
    PointType *PNearest = PNearestArray[ 0 ];

    t0 = GetCurrentTime();
    for ( std::size_t i = 0; i < 100; i++ ) {
      nodes_tree.SearchNearestPoint( allPoints, numsearch, PNearestArray, distancesArrayNearest );
    }
    t1 = GetCurrentTime();

    std::cout << t1 - t0 << "\t";
    double distance = 0;


    t0 = GetCurrentTime();
    for ( std::size_t i = 0; i < numsearch_nearest; i++ ) {
      PNearest = nodes_tree.SearchNearestPoint( point0, distance );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    std::cout << n << "\t" << *PNearestArray[ 0 ] << std::endl;
  };
}
