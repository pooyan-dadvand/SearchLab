/* -*- c++ -*- */
#pragma once

#include <iostream>
#include <vector>

#include "sphere_object.h"
#include "timer.h"

#include "configures/sphere_object_configure.h"
#include "containers.h"
#include "spatial_containers/spatial_containers.h"

namespace ObjectsOld {

  template < class SearchStructureType, class ObjectType, class IteratorType, class DistanceIterator >
    static void RunTests( char const Title[], IteratorType PBegin, IteratorType PEnd,
			  IteratorType results0, DistanceIterator distances0, std::size_t MaxResults,
			  ObjectType *allObjects, double radius0, std::size_t numsearch,
			  std::size_t bucket_size ) {
    static bool first_time = true;
    double t0, t1;
    std::size_t n;

    std::size_t npoints = PEnd - PBegin;
    std::vector< std::size_t > numresArray( numsearch );
    Entities::PtrObjectType objectToSearch = &allObjects[ 0 ];

    std::size_t numsearch_nearest = 100 * numsearch;

    t0 = GetCurrentTime();
    Kratos::BinsObjectDynamic< SphereObjectConfigure > nodes_tree( PBegin, PEnd );
    t1 = GetCurrentTime();

    // Header of timmings
    if ( first_time ) {
      std::cout << "SS\t\tGEN\tSIROMP\tSIRSER\tSNPOMP\tSNPSER\tNOFR\tNP" << std::endl;
      first_time = false;
    }
    std::cout << Title << "\t" << t1 - t0 << "\t";

    std::vector< std::vector< Entities::PtrObjectType > > objectResultsTmp;
#pragma omp parallel
    {
#pragma omp single
      {
	objectResultsTmp = std::vector< std::vector< Entities::PtrObjectType > >(
										 omp_get_num_threads(), std::vector< Entities::PtrObjectType >( npoints ) );
      }
    }

    t0 = GetCurrentTime();
#pragma omp parallel for
    for ( std::size_t i = 0; i < numsearch; i++ ) {
      std::vector< Entities::PtrObjectType >::iterator resultsItr =
        objectResultsTmp[ omp_get_thread_num() ].begin();
      DistanceIterator distanceItr = distances0;
      n = nodes_tree.SearchObjectsInRadius( objectToSearch, radius0, resultsItr, distanceItr,
					    MaxResults );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    t0 = GetCurrentTime();
    for ( std::size_t i = 0; i < numsearch; i++ ) {
      IteratorType resultsItr = results0;
      DistanceIterator distanceItr = distances0;
      n = nodes_tree.SearchObjectsInRadius( objectToSearch, radius0, resultsItr, distanceItr,
					    MaxResults );
    }
    t1 = GetCurrentTime();
    std::cout << t1 - t0 << "\t";

    ObjectType **PNearestArray = new ObjectType *[ numsearch ];
    std::vector< double > distancesArrayNearest( numsearch );
    ObjectType *PNearest = PNearestArray[ 0 ];

    // t0 = GetCurrentTime();
    // for (std::size_t i = 0; i < 100; i++) {
    // 	nodes_tree.SearchNearestPoint(allPoints, numsearch, PNearestArray, distancesArrayNearest);
    // }
    // t1 = GetCurrentTime();

    std::cout << "NaN"
	      << "\t";
    double distance = 0;


    // t0 = GetCurrentTime();
    // for (std::size_t i = 0; i < numsearch_nearest; i++) {
    // 	PNearest = nodes_tree.SearchNearestPoint(point0, distance);
    // }
    // t1 = GetCurrentTime();
    std::cout << "NaN"
	      << "\t";

    std::cout << n << "\t"
	      << "Does not apply" << std::endl;
  };

} // Namespace ObjectsOld
