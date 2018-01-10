#include <random>
#include "bins_cells_container.h"
#include "clock.h"

inline void GetSuggestedGridSize( size_t gridFull[ 3], const std::vector< Point> &points_vector) {
  gridFull[ 0] = 0;
  gridFull[ 1] = 0;
  gridFull[ 2] = 0;
  Crono clk;
  std::cout << "===== Calculate suggested grid size: ";
  BinsCellsContainer mCells_stats( points_vector.begin(), points_vector.end(), gridFull, false );
  float t = clk.fin();
  mCells_stats.PrintGridSize();
  gridFull[ 0 ] = mCells_stats.GetNumberOfCells( 0 );
  gridFull[ 1 ] = mCells_stats.GetNumberOfCells( 1 );
  gridFull[ 2 ] = mCells_stats.GetNumberOfCells( 2 );
  std::cout << "     time = " << t << "s." << std::endl;
}

inline void DoBinsStatistics( const std::vector< Point> &points_vector,
                              int discretization, bool select_random_points,
                              const size_t grid_sample[ 3], const char *filename) {
  std::vector< Point> points_sample;
  size_t npoints = points_vector.size();
  std::cout << "===== Occupancy statistics 1/" << discretization 
	    << " " << ( select_random_points ? "random" : "not random") 
	    << " points: =====" << std::endl;
  Crono clk_total;
  {
    Crono clk;
    if ( select_random_points ) {
      // select 1 out of discretization randomly, may be repeated points
      size_t npoints_sample = ( npoints + discretization - 1) / discretization;
      std::random_device rd;
      for ( size_t i = 0; i < npoints_sample; i++ ) {
	size_t idx = ( size_t )( rd() & 0x7ffffffff ) % npoints; // rd() is an unsigned int
	points_sample.push_back( points_vector[ idx ] );
      }
    } else {
      // select 1 out of discretization
      for ( size_t i = 0; i < npoints; i++ ) {
	if ( ( i % discretization ) == 0 )
	  points_sample.push_back( points_vector[ i ] );
      }
    }
    // float t = clk.fin();
    // std::cout << "     sampling time = " << t << "s." << std::endl;
  }

  {
    Crono clk;
    // size_t grid_sample[ 3 ] = { 0, 0, 0 }; // = { 55, 55, 55 };
    BinsCellsContainer cells_sample( points_sample.begin(), points_sample.end(), grid_sample );
    // float t = clk.fin();
    cells_sample.PrintStatistics();
    // std::cout << "     bins time = " << t << "s." << std::endl;
    cells_sample.PrintDensitiesInFile( filename);
  }

  float t = clk_total.fin();
  std::cout << "     statistics time = " << t << "s." << std::endl;
}
