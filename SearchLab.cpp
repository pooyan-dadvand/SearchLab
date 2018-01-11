// System includes
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include "win_missing.h"

#include <random>

std::size_t G_TotalNumberOfPoints = 0;
inline std::size_t GetTotalNumberOfPoints() { return G_TotalNumberOfPoints;}
inline void SetTotalNumberOfPoints( std::size_t num) { G_TotalNumberOfPoints = num;}

#ifdef USE_KRATOS
// Ugly fixes
#include <assert.h>

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

#include "clock.h"
#include "bins_statistics.h"

bool G_UsePointsBinsHash = false;
bool G_PrintBinsStatistics = false;
#ifdef NDEBUG
// i.e. compiling without assert, i.e. release mode
std::size_t G_NumberOfRepetitions = 1000000;
#else // NDEBUG
// compiling with asserts, i.e. debug mode
std::size_t G_NumberOfRepetitions = 1000;
#endif // NDEBUG

std::size_t G_GridSize[ 3] = { 0, 0, 0};

inline bool FileExists( const std::string &file_name) {
  // from https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
#ifndef _WIN32
  struct stat buffer;   
  int ret = stat( file_name.c_str(), &buffer ); // stat() fails on windows with files >
#else // _WIN32
  struct _stat64 buffer64;
  int ret = _stat64( file_name.c_str(), &buffer64 );
#endif // _WIN32
  return ( ret == 0); 
}

int RunPointSearchComparison( std::string Filename, double Radius ) {
  // Input data
  std::cout << std::setprecision( 4 ) << std::fixed;

  double t0 = GetCurrentTime();
#ifdef _WIN32
  std::ifstream input;
  input.open( Filename.c_str() );
  std::cout << "Gzip compressed streams are not supported in MS Windows\n"; // zstr has errors for bigger files ( may be > 4GB)
#else
  zstr::ifstream input( Filename.c_str() );
#endif

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

  // Point **points;
  Point point;
  // Object Interfaces ( not implemented at the moment )
  // SphereObject< 3 > object;
  // points = new Point *[ npoints ];
  std::vector< Point > points_vector;
  // Object Interfaces
  // - New Interface
  // TO BE FILLED
  // - Old Interface
#ifdef USE_KRATOS
  std::vector< Entities::PtrObjectType > objects( npoints );
#endif

  std::size_t pid;

  std::locale prev_loc = std::cout.getloc();
  std::cout.imbue( std::locale( "" ) ); // for thousand separators ...
  for ( std::size_t i = 0; i < npoints; i++ ) {
    input >> pid;
    input >> point;

    // if( !i || i == npoints - 1) {
    //   std::cout << "Point " << i + 1 << " = " << pid << " - " << point << std::endl;
    // }

    // Object Interfaces ( not implemented at the moment )
    // for ( std::size_t d = 0; d < 3; d++ ) {
    //   object[ d ] = point[ d ];
    // }
    // object.radius = 0.5 / ( double)npoints;

    // points[ i ] = new Point( point );
    // points[ i ]->id = pid;
    point.id = pid;
    points_vector.push_back( point);

    if ( ( ( i + 1) % 1000000 ) == 0 ) {
      std::cout << "     points read: " << i + 1 << " of " << npoints << std::endl;
    }

#ifdef USE_KRATOS
    // Object Interfaces ( not implemented at the moment )
    objects[ i ] = new SphereObject< 3 >( object );
    objects[ i ]->id = pid;
    objects[ i ]->radius = 0.5 / ( double)npoints;
#endif
  }
  std::cout.imbue( prev_loc ); // restore previous locale, i.e. without thousand separators
  double t1 = GetCurrentTime();
  std::cout << "Reading file = " << t1 - t0 << " sec." << std::endl;

  Point min_point( points_vector[ 0 ] );
  Point max_point( points_vector[ 0 ] );
  Point mid_point;
  SphereObject< 3 > mid_object;

  min_point.id = 0;
  max_point.id = 0;
  mid_point.id = 0;

  SetTotalNumberOfPoints( npoints);

  for ( std::size_t i = 0; i < npoints; i++ ) {
    for ( std::size_t j = 0; j < 3; j++ ) {
      if ( min_point[ j ] > ( points_vector[ i ] )[ j ] )
        min_point[ j ] = ( points_vector[ i ] )[ j ];
      if ( max_point[ j ] < ( points_vector[ i ] )[ j ] )
        max_point[ j ] = ( points_vector[ i ] )[ j ];
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

  std::size_t numrepetitions = G_NumberOfRepetitions;
  std::size_t numrepetitions_nearest = numrepetitions * 10;

  prev_loc = std::cout.getloc();
  std::cout.imbue( std::locale( "" ) ); // for thousand separators ...
  if ( G_PrintBinsStatistics) {
    std::cout << " min point : " << min_point << std::endl;
    std::cout << " max point : " << max_point << std::endl;
    std::cout << " search_point : " << search_point << std::endl;
    std::cout << " search radius : " << Radius << std::endl;
    std::cout << std::endl;
  }
  std::cout << " Number of Points : " << npoints << std::endl;
  std::cout << " Number of Repetitions : " << numrepetitions << std::endl;
  std::cout << std::endl;
  std::cout.imbue( prev_loc ); // restore previous locale, i.e. without thousand separators

  // Data Setup
  Point *allPoints = new Point[ numrepetitions ];
  SphereObject< 3 > *allSpheres = new SphereObject< 3 >[ numrepetitions ];

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

  bool do_guess_occupancy = true;
  char filename_pnts[ 10240];
  *filename_pnts = '\0';
  if ( G_PrintBinsStatistics && do_guess_occupancy ) {
    Crono clk;
    size_t gridFull[ 3 ] = { 0, 0, 0 };
    GetSuggestedGridSize( gridFull, points_vector);

    char filename_scr[ 10240];
    char filename_png[ 10240];
    strcpy( filename_pnts, Filename.c_str());
    char *ext = strrchr( filename_pnts, '.');
    if ( !ext)
      ext = &filename_pnts[ strlen( filename_pnts)];
    strcpy( ext, "_gnuplot.dat");
    strcpy( filename_scr, Filename.c_str());
    ext = strrchr( filename_scr, '.');
    if ( !ext)
      ext = &filename_scr[ strlen( filename_scr)];
    strcpy( ext, "_gnuplot.cmd");
    strcpy( filename_png, Filename.c_str());
    ext = strrchr( filename_png, '.');
    if ( !ext)
      ext = &filename_png[ strlen( filename_png)];
    strcpy( ext, "_gnuplot.png");

    FILE *fo = fopen( filename_scr, "w");
    fprintf( fo, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5   # --- blue\n");
    fprintf( fo, "set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 1.5   # --- red\n");
    fprintf( fo, "set style line 3 lc rgb '#18dd1f' lt 1 lw 2 pt 5 ps 1.5   # --- green\n");
    fprintf( fo, "set style line 4 lc rgb '#dddd1f' lt 1 lw 2 pt 5 ps 1.5   # --- yellow\n");
    fprintf( fo, "set logscale xy 10\n");
    fprintf( fo, "set title '%s ( normalized )'\n", Filename.c_str());
    fprintf( fo, "plot '%s' using 1:3 index 0 with linespoints ls 1 title '1/1000', ", filename_pnts);
    fprintf( fo, " '' using 1:3 index 1 with linespoints ls 2 title '1/100', ");
    fprintf( fo, " '' using 1:3 index 2 with linespoints ls 3 title '1/10', ");
    fprintf( fo, " '' using 1:3 index 3 with linespoints ls 4 title 'full' \n");
    fprintf( fo, "set term png size 1280,720\n");
    fprintf( fo, "set output '%s'\n", filename_png);
    fprintf( fo, "replot\n");
    fprintf( fo, "set term x11\n");
    fclose( fo);

    fo = fopen( filename_pnts, "w");
    fprintf( fo, "# cell densities information for '%s'\n\n", Filename.c_str());
    fclose( fo);

    int discretization = 1000; // get_ceil_prime( 1000);
    bool select_random_points = false;
    gridFull[ 0 ] = 0;
    gridFull[ 1 ] = 0;
    gridFull[ 2 ] = 0;
    // DoBinsStatistics( points_vector, discretization,  select_random_points, gridFull );
    select_random_points = true;
    DoBinsStatistics( points_vector, discretization,  select_random_points, gridFull, filename_pnts );
    discretization = 100; //get_ceil_prime( 100);
    // select_random_points = false;
    // DoBinsStatistics( points_vector, discretization,  select_random_points, gridFull );
    select_random_points = true;
    DoBinsStatistics( points_vector, discretization,  select_random_points, gridFull, filename_pnts );
    discretization = 10;//get_ceil_prime( 10);
    // select_random_points = false;
    // DoBinsStatistics( points_vector, discretization,  select_random_points, gridFull );
    select_random_points = true;
    DoBinsStatistics( points_vector, discretization,  select_random_points, gridFull, filename_pnts );
    
    float t = clk.fin();
    std::cout << "     statistics time = " << t << "s." << std::endl;
    std::cout << "Created '" << filename_pnts << "' and '" << filename_scr << "'" << std::endl;
    std::cout << "To generate '" << filename_png << "' try:" << std::endl;
    std::cout << "     gnuplot " << filename_scr << std::endl;
  }

  // Point Interfaces
  // - New Interface
  if ( !G_UsePointsBinsHash) {
    std::cout << "Using PointsBins" << std::endl;
    PointsNew::RunTests< PointsBins< Point > >( "PointBins", points_vector, search_point, Radius,
						numrepetitions, numrepetitions_nearest,
						G_GridSize, G_PrintBinsStatistics, filename_pnts);
  } else {
    std::cout << "Using PointsHash" << std::endl;
    PointsNew::RunTests< PointsBinsHash< Point > >( "PointBinsHash", points_vector, search_point, Radius,
						    numrepetitions, numrepetitions_nearest,
						    G_GridSize, G_PrintBinsStatistics, filename_pnts );
  }

// - Old Interface
#ifdef USE_KRATOS
  PointsOld::RunTests< Containers::BinsStaticType >( "StaticBins", points, points + npoints,
                                                     p_results, resultDistances.begin(),
                                                     max_results, allPoints, Radius, numrepetitions, 1 );
  PointsOld::RunTests< Containers::BinsDynamicType >(
      "DynamicBins", points, points + npoints, p_results, resultDistances.begin(), max_results,
      allPoints, Radius, numrepetitions, 1 );
  PointsOld::RunTests< Containers::OctreeType >( "OctTree\t", points, points + npoints, p_results,
                                                 resultDistances.begin(), max_results, allPoints,
                                                 Radius, numrepetitions, 10 );
  PointsOld::RunTests< Containers::BinsDynamicOctreeType >(
      "OcTreeDynamic\t", points, points + npoints, p_results, resultDistances.begin(), max_results,
      allPoints, Radius, numrepetitions, 10 );
#endif

// Object Interfaces
// - New Interface
// TO BE FILLED

// - Old Interface
#ifdef USE_KRATOS
  ObjectsOld::RunTests< Containers::BinsObjectStaticType >(
      "StaticObjects", objects.begin(), objects.end(), objectResults.begin(),
      resultDistances.begin(), max_results, allSpheres, Radius, numrepetitions, 1 );
  ObjectsOld::RunTests< Containers::BinsObjectDynamicType >(
      "DynamicObjects", objects.begin(), objects.end(), objectResults.begin(),
      resultDistances.begin(), max_results, allSpheres, Radius, numrepetitions, 1 );
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

void PrintUsage( const std::string &progname, const std::string &filename, double radius) {
  std::cout << "Usage: " << progname << " [ -type bin/hash] [ -num num_repetitions] [ -grid XxYxZ] [ -statistics] [ filename [ radius ] ]" << std::endl;
  std::cout << "       filename default value = " << filename << std::endl;
  std::cout << "       radius default value   = " << radius << std::endl;
  std::cout << "       type default value = " << ( G_UsePointsBinsHash ? "hash" : "bin") << std::endl;
  std::cout << "       num_repetitions default value = " << G_NumberOfRepetitions << std::endl;
  std::cout << "       grid size (XxYxZ) default depends on number of points" << std::endl;
  std::cout << "       statistics default value is off" << std::endl;
}

int main( int argc, char *argv[] ) {

  // testParallelCoherentHash();
  // return 0;

    
  // Default filename
  std::string filename_default = "../cases/genericCube2x2x2.500000.pts";
  // filename = "../cases/genericCube100x100x100.5051.pts";
  // filename = "../cases/randomCube2000000.pts";

  // Default radius
  double radius_default = 0.0102;

  // used filename and radius
  std::string filename = filename_default;
  double radius = radius_default;

  if ( argc > 1 ) {
    bool has_filename = false;
    for ( int idx_arg = 1; idx_arg < argc; idx_arg++) {
      if ( argv[ idx_arg][ 0] == '-') {
	if ( !strncasecmp( argv[ idx_arg ], "-h", 2 ) || !strncasecmp( argv[ idx_arg ], "--h", 3 ) ) {
	  PrintUsage( argv[ 0], filename_default, radius_default);
	  return 0;
	} else if ( !strncasecmp( argv[ idx_arg ], "-t", 2 ) || !strncasecmp( argv[ idx_arg ], "--t", 3 ) ) {
	  if ( idx_arg + 1 < argc) {
	    idx_arg++;
	    if ( !strncasecmp( argv[ idx_arg ], "b", 1 ) || !strcasecmp( argv[ idx_arg ], "bin") ) {
	      G_UsePointsBinsHash = false;
	    } else if ( !strncasecmp( argv[ idx_arg ], "h", 1 ) || !strcasecmp( argv[ idx_arg ], "hash" ) ) {
	      G_UsePointsBinsHash = true;
	    } else {
	      std::cout << "error: unknown type: " << argv[ idx_arg] << std::endl;;
	      PrintUsage( argv[ 0], filename_default, radius_default);
	      return 0;
	    }
	  } else {
	    std::cout << "error: missing type" << std::endl;;
	    PrintUsage( argv[ 0], filename_default, radius_default);
	    return 0;
	  }
	} else if ( !strncasecmp( argv[ idx_arg ], "-n", 2 ) || !strncasecmp( argv[ idx_arg ], "--n", 3 ) ) {
	  if ( idx_arg + 1 < argc) {
	    idx_arg++;
	    size_t numrepetitions = 0;
	    if ( sscanf( argv[ idx_arg], "%zu", &numrepetitions) == 1) {
	      G_NumberOfRepetitions = numrepetitions;
	    } else {
	      std::cout << "error: unknown number of searches: " << argv[ idx_arg] << std::endl;;
	      PrintUsage( argv[ 0], filename_default, radius_default);
	      return 0;
	    }
	  } else {
	    std::cout << "error: missing number of searches" << std::endl;;
	    PrintUsage( argv[ 0], filename_default, radius_default);
	    return 0;
	  }
	} else if ( !strncasecmp( argv[ idx_arg ], "-g", 2 ) || !strncasecmp( argv[ idx_arg ], "--g", 3 ) ) {
	  if ( idx_arg + 1 < argc) {
	    idx_arg++;
	    size_t grid_size_x = 0;
	    size_t grid_size_y = 0;
	    size_t grid_size_z = 0;
	    if ( sscanf( argv[ idx_arg], "%zux%zux%zu", &grid_size_x, &grid_size_y, &grid_size_z) == 3) {
	      G_GridSize[ 0] = grid_size_x;
	      G_GridSize[ 1] = grid_size_y;
	      G_GridSize[ 2] = grid_size_z;
	    } else {
	      std::cout << "error: unknown grid sizes (XxYxZ): " << argv[ idx_arg] << std::endl;;
	      PrintUsage( argv[ 0], filename_default, radius_default);
	      return 0;
	    }
	  } else {
	    std::cout << "error: missing grid sizes (XxYxZ)" << std::endl;;
	    PrintUsage( argv[ 0], filename_default, radius_default);
	    return 0;
	  }
	} else if ( !strncasecmp( argv[ idx_arg ], "-s", 2 ) || !strncasecmp( argv[ idx_arg ], "--s", 3 ) ) {
	  G_PrintBinsStatistics  = true;
	}
      } else {
	if ( !has_filename) {
	  filename = argv[ idx_arg ];
	  has_filename = true;
	}
	radius = atof( argv[ idx_arg ] ) / 1000000;
	break; // the last argument is the radius if present
      }
    }
  }

  if ( !FileExists( filename)) {
    std::cout << "File '" << filename << "' can not be accessed." << std::endl;
    PrintUsage( argv[ 0], filename_default, radius_default);
    return 0;
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

// #include "clock.cc"
