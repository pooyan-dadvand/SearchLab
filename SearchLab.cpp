// System includes
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <win_missing.h>

#ifdef USE_KRATOS
// Ugly fixes
#include <assert.h>
#define KRATOS_ERROR std::cout

#include "interfaces/points_old_interface.h"
#include "interfaces/objects_old_interface.h"
#endif

// Support for compressed streams
#include "zstr.hpp"

// Containers
#include "containers.h"

// Interfaces
#include "interfaces/points_new_interface.h"

#include "parallel_coherent_hash.h"

int RunPointSearchComparison( std::string Filename, double Radius ) {

  // Input data
  std::cout << std::setprecision( 4 ) << std::fixed;

  Point **points;
  Point point;
  SphereObject< 3 > object;

  double t0 = GetCurrentTime();
  // std::ifstream input;
  // input.open( Filename.c_str() );
  zstr::ifstream input( Filename.c_str() );

  if ( !input ) {
    std::cout << "Cannot open data file" << std::endl;
    return 1;
  }

  std::cout << "Comparison for " << Filename << std::endl;

  std::size_t npoints;

  input >> npoints;

  if ( npoints <= 0) {
    std::cout << "ERROR: no points information read, nothing to do!!!" << std::endl;
    return 1;
  }
  
  points = new Point *[ npoints ];
  std::vector< Entities::PtrObjectType > objects( npoints );

  std::size_t pid;

  for ( std::size_t i = 0; i < npoints; i++ ) {
    input >> pid;
    input >> point;

    // if( !i || i == npoints - 1) {
    //   std::cout << "Point " << i + 1 << " = " << pid << " - " << point << std::endl;
    // }

    for ( std::size_t d = 0; d < 3; d++ ) {
      object[ d ] = point[ d ];
    }
    object.radius = 0.5 / ( double)npoints;

    points[ i ] = new Point( point );
    points[ i ]->id = pid;

    objects[ i ] = new SphereObject< 3 >( object );
    objects[ i ]->id = pid;
    objects[ i ]->radius = 0.5 / ( double)npoints;
  }
  double t1 = GetCurrentTime();
  std::cout << "Reading file = " << t1 - t0 << " sec." << std::endl;

  Point min_point( *points[ 0 ] );
  Point max_point( *points[ 0 ] );
  Point mid_point;
  SphereObject< 3 > mid_object;

  min_point.id = 0;
  max_point.id = 0;
  mid_point.id = 0;

  for ( std::size_t i = 0; i < npoints; i++ ) {
    for ( std::size_t j = 0; j < 3; j++ ) {
      if ( min_point[ j ] > ( *points[ i ] )[ j ] )
        min_point[ j ] = ( *points[ i ] )[ j ];
      if ( max_point[ j ] < ( *points[ i ] )[ j ] )
        max_point[ j ] = ( *points[ i ] )[ j ];
    }
  }

  for ( std::size_t i = 0; i < Dim; i++ ) {
    mid_point.coord[ i ] = ( max_point[ i ] + min_point[ i ] ) / 2.00;
    mid_object.coord[ i ] = ( max_point[ i ] + min_point[ i ] ) / 2.00;
  }

  mid_object.radius = 0.5 / ( double)npoints;

  // Output data Info
  Point &search_point = mid_point;
  SphereObject< 3 > &search_object = mid_object;

  std::size_t numsearch = 1000000;
  std::size_t numsearch_nearest = numsearch * 10;

  std::cout << " min point : " << min_point << std::endl;
  std::cout << " max point : " << max_point << std::endl;
  std::cout << " search_point : " << search_point << std::endl;
  std::cout << " search radius : " << Radius << std::endl;
  std::cout << std::endl;

  std::cout << " Number of Points : " << npoints << std::endl;
  std::cout << " Number of Repetitions : " << numsearch << std::endl;
  std::cout << std::endl;

  // Data Setup
  Point *allPoints = new Point[ numsearch ];
  SphereObject< 3 > *allSpheres = new SphereObject< 3 >[ numsearch ];

  std::size_t max_results = npoints;
  for ( std::size_t i = 0; i < 1; i++ ) {
    allPoints[ i ] = search_point;
    allSpheres[ i ] = search_object;
  }

  // Prepare the search point, search radius and resut arrays

  std::vector< Entities::PtrObjectType > objectResults( max_results );
  std::vector< double > resultDistances( max_results );

#ifdef USE_KRATOS
  double *distances = new double[ npoints ];
  Entities::PointIterator p_results = new Entities::PtrPointType[ max_results ];
#endif // USE_KRATOS

  // Point-Based Search Structures
  std::vector< Point > points_vector;
  for ( std::size_t i = 0; i < npoints; i++ ) {
    points_vector.push_back( *( points[ i ] ) );
  }

  // Point Interfaces
  // - New Interface
  PointsNew::RunTests< PointsBins< Point > >( "PointBins", points_vector, search_point, Radius,
                                              numsearch, numsearch_nearest );
  PointsNew::RunTests< PointsBinsHash< Point > >( "PointBins", points_vector, search_point, Radius,
						  numsearch, numsearch_nearest );

// - Old Interface
#ifdef USE_KRATOS
  PointsOld::RunTests< Containers::BinsStaticType >( "StaticBins", points, points + npoints,
                                                     p_results, resultDistances.begin(),
                                                     max_results, allPoints, Radius, numsearch, 1 );
  PointsOld::RunTests< Containers::BinsDynamicType >(
      "DynamicBins", points, points + npoints, p_results, resultDistances.begin(), max_results,
      allPoints, Radius, numsearch, 1 );
  PointsOld::RunTests< Containers::OctreeType >( "OctTree\t", points, points + npoints, p_results,
                                                 resultDistances.begin(), max_results, allPoints,
                                                 Radius, numsearch, 10 );
  PointsOld::RunTests< Containers::BinsDynamicOctreeType >(
      "OcTreeDynamic\t", points, points + npoints, p_results, resultDistances.begin(), max_results,
      allPoints, Radius, numsearch, 10 );
#endif

// Object Interfaces
// - New Interface
// TO BE FILLED

// - Old Interface
#ifdef USE_KRATOS
  ObjectsOld::RunTests< Containers::BinsObjectStaticType >(
      "StaticObjects", objects.begin(), objects.end(), objectResults.begin(),
      resultDistances.begin(), max_results, allSpheres, Radius, numsearch, 1 );
  ObjectsOld::RunTests< Containers::BinsObjectDynamicType >(
      "DynamicObjects", objects.begin(), objects.end(), objectResults.begin(),
      resultDistances.begin(), max_results, allSpheres, Radius, numsearch, 1 );
#endif
  // RunTestsOldInterface<BinsObjectDynamicType>

  return 0;
}

void testParallelCoherentHash() {
  ParallelCoherentHash< int, int> pch;
  const int size_hash = 1000;
  pch.resize( size_hash, -1);
  for ( int i = 0; i < size_hash; i++) {
    // pch[ i] = i + 1;
    int key = i - 5;
    pch.getDataRef( key) = i + 1;
  }
  for ( int i = 0; i < size_hash; i++) {
    int key = i - 5;
    pch.getDataRef( key)++;
    // printf( "id[ %4d] = %d\n", key, pch.getData( key));
  }
  pch.PrintStatistics();
}

int main( int arg, char *argv[] ) {

  // testParallelCoherentHash();
  // return 0;

    
  std::string filename;

  // Default filename
  filename = "../cases/genericCube2x2x2.500000.pts";
  // filename = "../cases/genericCube100x100x100.5051.pts";
  // filename = "../cases/randomCube2000000.pts";

  // Default radius
  double radius = 0.0102;

  if ( arg > 1 ) {
    if ( !strncasecmp( argv[ 1 ], "-h", 2 ) || !strncasecmp( argv[ 1 ], "--h", 2 ) ) {
      std::cout << "Usage: " << argv[ 0 ] << " [ filename [ radius ] ]" << std::endl;
      std::cout << "       filename default value = " << filename << std::endl;
      std::cout << "       radius default value   = " << radius << std::endl;
      return 0;
    }
    filename = argv[ 1 ];
    if ( arg == 3 ) {
      radius = atof( argv[ 2 ] ) / 1000000;
    }
  }

  RunPointSearchComparison( filename, radius );

  // filename = "../cases/offsetCube79x79x79.1603.pts";
  // RunPointSearchComparison(filename, radius);
  // filename = "../cases/clusterCube6x6x6X4913.490.pts";
  // RunPointSearchComparison(filename, radius);
  // filename = "../cases/line100000.5.pts";
  // RunPointSearchComparison(filename, radius);

  return 0;
}
