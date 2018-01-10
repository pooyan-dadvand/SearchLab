/* -*- c++ -*- */
#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include <iostream>
#include <locale>

#include "bins_cells_container_hash.h"
#include "spatial_search_result.h"

template<typename T>
T abs_diff( T a, T b ) {
  return a > b ? a - b : b - a;
}

template < typename TObjectType > class PointsBinsHash {
  static constexpr int Dimension = 3;

public:
  using InternalPointType = std::array< double, Dimension >;
  using ObjectType = TObjectType;
  using PointerType = TObjectType *;
  using ResultType = SpatialSearchResult< TObjectType >;

  template < typename TIteratorType >
  PointsBinsHash( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd,
		  const std::size_t GridSize[ 3] )
    : mCells( PointsBegin, PointsEnd, GridSize ) {
    mNumberOfPoints = std::distance( PointsBegin, PointsEnd );

    if ( mNumberOfPoints == 0 ) {
      mpPoints = nullptr;
      return;
    }
    mpPoints = new PointerType[ mNumberOfPoints ];
    for ( std::size_t i = 0; i < mNumberOfPoints; i++ )
      mpPoints[ i ] = nullptr;

    AssignPointsToCells( PointsBegin, PointsEnd );
  }

  void SearchInRadius( TObjectType const &ThePoint, double Radius,
                       std::vector< ResultType > &rResults ) {
    InternalPointType min_point;
    std::array< std::size_t, Dimension > length;
    const double radius2e = std::pow( Radius + std::numeric_limits< double >::epsilon(), 2.0 );

    for ( int i = 0; i < Dimension; i++ ) {
      min_point[ i ] = ThePoint[ i ] - Radius;
      length[ i ] = mCells.CalculatePosition( ThePoint[ i ] + Radius, i ) -
                    mCells.CalculatePosition( ThePoint[ i ] - Radius, i ) + 1;
    }
    auto min_cell = mCells.CalculateCellIndex( min_point );

    for ( std::size_t i_z = 0; i_z < length[ 2 ]; i_z++ ) {
      auto y_position =
          min_cell + i_z * mCells.GetNumberOfCells( 0 ) * mCells.GetNumberOfCells( 1 );
      for ( std::size_t i_y = 0; i_y < length[ 1 ]; i_y++ ) {
	std::size_t offset, end_offset;
	bool found = mCells.GetCellStoredOffsets( y_position, offset, end_offset);
	if ( !found) continue;

	TObjectType **p_point = mpPoints + offset;
	for ( ; offset < end_offset; offset++ ) {
	  if ( Distance2( **p_point, ThePoint ) <= radius2e ) {
	    rResults.push_back( ResultType( *p_point ) );
	  }
	  p_point++;
	}
	y_position += mCells.GetNumberOfCells( 0 );
      }
    }
  }

  ResultType SearchNearest( TObjectType const &ThePoint ) {
    return this->NewSearchNearest( ThePoint );

    // auto cell_index = mCells.CalculateCellIndex(ThePoint);
    ResultType current_result;

    if ( mNumberOfPoints == 0 )
      return current_result;

    current_result.SetDistance2( std::numeric_limits< double >::max() );
    double radius = std::max( mCells.GetCellSize( 0 ), mCells.GetCellSize( 1 ) );
    radius = std::max( radius, mCells.GetCellSize( 2 ) ) * .5;

    while ( !current_result.IsObjectFound() ) {
      InternalPointType min_point;
      std::array< std::size_t, Dimension > length;
      for ( int i = 0; i < Dimension; i++ ) {
        min_point[ i ] = ThePoint[ i ] - radius;
        length[ i ] = mCells.CalculatePosition( ThePoint[ i ] + radius, i ) -
                      mCells.CalculatePosition( ThePoint[ i ] - radius, i ) + 1;
      }
      auto min_cell = mCells.CalculateCellIndex( min_point );

      for ( std::size_t i_z = 0; i_z < length[ 2 ]; i_z++ ) {
        auto y_position =
            min_cell + i_z * mCells.GetNumberOfCells( 0 ) * mCells.GetNumberOfCells( 1 );
        for ( std::size_t i_y = 0; i_y < length[ 1 ]; i_y++ ) {
	  for ( std::size_t i_x = 0; i_x < length[ 2 ]; i_x++) {
	    auto x_position = y_position + i_x;
	    std::size_t offset, end_offset;
	    bool found = mCells.GetCellStoredOffsets( x_position, offset, end_offset);
	    if ( !found) continue;

	    TObjectType **p_point = mpPoints + offset;//mCells.GetCellBeginIndex( y_position );
	    for ( ; offset < end_offset; offset++ ) {
	      double distance2 = Distance2( **p_point, ThePoint );
	      if ( distance2 < current_result.GetDistance2() ) {
		current_result.Set( *p_point );
		current_result.SetDistance2( distance2 );
	      }
	      p_point++;
	    }
          } // for i_x
          y_position += mCells.GetNumberOfCells( 0 );
        } // for i_y
      } // for i_z
      radius *= 2.00;
    }

    return current_result;
  }

  ResultType SearchNearestWithinRadius( TObjectType const &ThePoint, double limit_radius ) {
    // auto cell_index = mCells.CalculateCellIndex(ThePoint);
    ResultType current_result;

    if ( mNumberOfPoints == 0 )
      return current_result;

    current_result.SetDistance2( std::numeric_limits< double >::max() );
    double radius = std::max( mCells.GetCellSize( 0 ), mCells.GetCellSize( 1 ) );
    radius = std::max( radius, mCells.GetCellSize( 2 ) ) * .5;

    while ( !current_result.IsObjectFound() && ( radius < limit_radius)) {
      InternalPointType min_point;
      std::array< std::size_t, Dimension > length;
      for ( int i = 0; i < Dimension; i++ ) {
        min_point[ i ] = ThePoint[ i ] - radius;
        length[ i ] = mCells.CalculatePosition( ThePoint[ i ] + radius, i ) -
                      mCells.CalculatePosition( ThePoint[ i ] - radius, i ) + 1;
      }
      auto min_cell = mCells.CalculateCellIndex( min_point );

      for ( std::size_t i_z = 0; i_z < length[ 2 ]; i_z++ ) {
        auto y_position =
            min_cell + i_z * mCells.GetNumberOfCells( 0 ) * mCells.GetNumberOfCells( 1 );
        for ( std::size_t i_y = 0; i_y < length[ 1 ]; i_y++ ) {
	  for ( std::size_t i_x = 0; i_x < length[ 2 ]; i_x++) {
	    auto x_position = y_position + i_x;
	    std::size_t offset, end_offset;
	    bool found = mCells.GetCellStoredOffsets( x_position, offset, end_offset);
	    if ( !found) continue;

	    TObjectType **p_point = mpPoints + offset;//mCells.GetCellBeginIndex( y_position );
	    for ( ; offset < end_offset; offset++ ) {
	      double distance2 = Distance2( **p_point, ThePoint );
	      if ( distance2 < current_result.GetDistance2() ) {
		current_result.Set( *p_point );
		current_result.SetDistance2( distance2 );
	      }
	      p_point++;
	    }
          } // for i_x
          y_position += mCells.GetNumberOfCells( 0 );
        } // for i_y
      } // for i_z
      radius *= 2.00;
    }
    return current_result;
  }
 

  ResultType NewSearchNearest( TObjectType const &ThePoint ) {
    ResultType current_result;
    if ( mNumberOfPoints == 0 )
      return current_result;

    current_result.SetDistance2( std::numeric_limits< double >::max() );
    double cell_diagonal = sqrt( mCells.GetCellSize( 0 ) * mCells.GetCellSize( 0 ) + 
				 mCells.GetCellSize( 1 ) * mCells.GetCellSize( 1 ) + 
				 mCells.GetCellSize( 2 ) * mCells.GetCellSize( 2 ));
    double limit_radius = 2.0 * cell_diagonal; // look into ThePoint cell and surrounding cells of the 1st shell

    current_result = this->SearchNearestWithinRadius( ThePoint, limit_radius);
    if ( current_result.IsObjectFound()) {
      return current_result;
    }

    // nearest point not found let's do a search over occupied cells
    const std::vector< std::size_t> &lstUsedCells = mCells.GetListUsedCellIndices();

    // first round to get the closest cell centre to ThePoint
    bool use_centre = true;
    double distance2_centre_nearest_cell = std::numeric_limits< double >::max();
    long long distance2_idx_nearest_cell = std::numeric_limits< long long >::max();
    if ( use_centre ) {
      for ( auto it_idx = lstUsedCells.begin(); it_idx < lstUsedCells.end(); it_idx++ ) {
        std::size_t idx = *it_idx;
        InternalPointType cell_centre = mCells.CalculateCentreOfCell( idx );
        double distance2 = Distance2( Point( cell_centre[ 0 ], cell_centre[ 1 ], cell_centre[ 2 ] ), ThePoint );
        if ( distance2 < distance2_centre_nearest_cell ) {
          distance2_centre_nearest_cell = distance2;
        }
      }
    } else {
      std::size_t idx_search_cell = mCells.CalculateCellIndex( ThePoint );
      std::size_t idx_search_x, idx_search_y, idx_search_z;
      std::size_t idx_cell_min = 0;
      mCells.GetCellVectorIndices( idx_search_cell, idx_search_x, idx_search_y, idx_search_z );
      for ( auto it_idx = lstUsedCells.begin(); it_idx < lstUsedCells.end(); it_idx++ ) {
        std::size_t idx = *it_idx;
        std::size_t idx_x, idx_y, idx_z;
        mCells.GetCellVectorIndices( idx, idx_x, idx_y, idx_z );
        long long diff_x = abs( ( long long )idx_x - ( long long )idx_search_x ); // abs_diff( idx_x, idx_search_x );
        long long diff_y = abs( ( long long )idx_y - ( long long )idx_search_y ); //abs_diff( idx_y, idx_search_y );
        long long diff_z = abs( ( long long )idx_z - ( long long )idx_search_z ); //abs_diff( idx_z, idx_search_z );
        long long dist2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
        if ( dist2 < distance2_idx_nearest_cell ) {
          distance2_idx_nearest_cell = dist2;
          idx_cell_min = idx;
        }
      }
      InternalPointType cell_centre = mCells.CalculateCentreOfCell( idx_cell_min );
      distance2_centre_nearest_cell = Distance2( Point( cell_centre[ 0 ], cell_centre[ 1 ], cell_centre[ 2 ] ), ThePoint );
    }
    
    if ( distance2_centre_nearest_cell == std::numeric_limits< double >::max()) // == no points in bin !!!
      return current_result;

    // now get the used cells within this point + radius
    if ( use_centre ) {
      double search_radius = distance2_centre_nearest_cell + cell_diagonal * 1.5;
      for ( auto it_idx = lstUsedCells.begin(); it_idx < lstUsedCells.end(); it_idx++ ) {
        std::size_t idx = *it_idx;
        InternalPointType cell_centre = mCells.CalculateCentreOfCell( idx );
        double distance2 = Distance2( Point( cell_centre[ 0 ], cell_centre[ 1 ], cell_centre[ 2 ] ), ThePoint );
        if ( distance2 < search_radius ) {
          // look into the points of the cell
          this->SearchNearestInCell( idx, ThePoint, current_result );
        }
      }
    } else {
      long long search_radius = distance2_idx_nearest_cell + 3; // i.e. diagonal of idx = sqrt( 1*1 + 1*1 + 1*1), so idx*idx = 3
      std::size_t idx_search_cell = mCells.CalculateCellIndex( ThePoint );
      std::size_t idx_search_x, idx_search_y, idx_search_z;
      mCells.GetCellVectorIndices( idx_search_cell, idx_search_x, idx_search_y, idx_search_z );
      for ( auto it_idx = lstUsedCells.begin(); it_idx < lstUsedCells.end(); it_idx++ ) {
        std::size_t idx = *it_idx;
        std::size_t idx_x, idx_y, idx_z;
        mCells.GetCellVectorIndices( idx, idx_x, idx_y, idx_z );
        long long diff_x = abs( ( long long )idx_x - ( long long )idx_search_x ); // abs_diff( idx_x, idx_search_x );
        long long diff_y = abs( ( long long )idx_y - ( long long )idx_search_y ); //abs_diff( idx_y, idx_search_y );
        long long diff_z = abs( ( long long )idx_z - ( long long )idx_search_z ); //abs_diff( idx_z, idx_search_z );
        long long dist2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
        if ( dist2 < search_radius ) {
          // look into the points of the cell
          this->SearchNearestInCell( idx, ThePoint, current_result );
        }
      }
    }

    return current_result;
  }

  void PrintStatistics() const {
    // Bins statistics
    mCells.PrintStatistics();
  }

private:
  std::size_t mNumberOfPoints;
  BinsCellsContainerHash mCells;
  TObjectType **mpPoints;

  template < typename TIteratorType >
  void AssignPointsToCells( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd ) {
    for ( auto i_point = PointsBegin; i_point != PointsEnd; i_point++ ) {
      auto index = mCells.CalculateCellIndex( *i_point );

      std::size_t offset_begin, offset_end;
      bool found = mCells.GetCellStoredOffsets( index, offset_begin, offset_end);
      if ( !found) continue;
      for ( std::size_t offset = offset_begin; offset < offset_end; offset++) {
        if ( mpPoints[ offset ] == nullptr ) {
          mpPoints[ offset ] = &( *i_point );
          break;
        }
      }
    }
  }

  void SearchNearestInCell( std::size_t CellIndex, TObjectType const &ThePoint,
                            ResultType &rCurrentResult ) {
    std::size_t offset_begin, offset_end;
    bool found = mCells.GetCellStoredOffsets( CellIndex, offset_begin, offset_end);
    if ( !found)
      return;
    for ( std::size_t offset = offset_begin; offset < offset_end; offset++ ) {
      TObjectType *p_point = mpPoints[ offset ];
      auto distance2 = Distance2( *p_point, ThePoint );
      if ( distance2 <= rCurrentResult.GetDistance2() ) {
        rCurrentResult.Set( p_point );
        rCurrentResult.SetDistance2( distance2 );
      }
    }
  }

  double Distance2( TObjectType const &FirstPoint, TObjectType const &SecondPoint ) {
    double result = double();
    for ( int i = 0; i < Dimension; i++ ) {
      auto distance_i = FirstPoint[ i ] - SecondPoint[ i ];
      result += distance_i * distance_i;
    }
    return result;
  }
};
