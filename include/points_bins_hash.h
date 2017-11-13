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

template < typename TObjectType > class PointsBinsHash {
  static constexpr int Dimension = 3;

public:
  using InternalPointType = std::array< double, Dimension >;
  using ObjectType = TObjectType;
  using PointerType = TObjectType *;
  using ResultType = SpatialSearchResult< TObjectType >;

  template < typename TIteratorType >
  PointsBinsHash( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd )
    : mCells( PointsBegin, PointsEnd ) {
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
	std::size_t offset = mCells.GetCellBeginIndex( y_position );
	// const std::size_t end_offset = mCells.GetCellBeginIndex( y_position + length[ 0 ] );
	const std::size_t end_offset = mCells.GetCellEndIndex( y_position);
	TObjectType **p_point = mpPoints + offset;//mCells.GetCellBeginIndex( y_position );
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
      // const double radius2 = radius * radius;
      // std::cout << " ... index search point = " << mCells.CalculateCellIndex( ThePoint) << std::endl;
      // std::cout << " ... index search point( 0) = " << mCells.CalculatePosition( ThePoint[ 0 ], 0) << std::endl;
      // std::cout << " ... index search point( 1) = " << mCells.CalculatePosition( ThePoint[ 1 ], 1) << std::endl;
      // std::cout << " ... index search point( 2) = " << mCells.CalculatePosition( ThePoint[ 2 ], 2) << std::endl;
      for ( int i = 0; i < Dimension; i++ ) {
        min_point[ i ] = ThePoint[ i ] - radius;
        length[ i ] = mCells.CalculatePosition( ThePoint[ i ] + radius, i ) -
                      mCells.CalculatePosition( ThePoint[ i ] - radius, i ) + 1;
      }
      // std::cout << " ... length = " << length[ 0] << " " << length[ 1] << " " << length[ 2] << std::endl;
      auto min_cell = mCells.CalculateCellIndex( min_point );

      for ( std::size_t i_z = 0; i_z < length[ 2 ]; i_z++ ) {
        auto y_position =
            min_cell + i_z * mCells.GetNumberOfCells( 0 ) * mCells.GetNumberOfCells( 1 );
        for ( std::size_t i_y = 0; i_y < length[ 1 ]; i_y++ ) {
	  for ( std::size_t i_x = 0; i_x < length[ 2 ]; i_x++) {
	    auto x_position = y_position + i_x;
	    // std::cout << " ... testing cell = " << x_position << std::endl;
	    std::size_t offset = mCells.GetCellBeginIndex( x_position );
	    const std::size_t end_offset = mCells.GetCellEndIndex( x_position); // mCells.GetCellBeginIndex( y_position + length[ 0 ] );
	    TObjectType **p_point = mpPoints + offset;//mCells.GetCellBeginIndex( y_position );
	    for ( ; offset < end_offset; offset++ ) {
	      double distance_2 = Distance2( **p_point, ThePoint );
	      // std::cout << " ... looking into point # " << ( *p_point)->id << std::endl;
	      if ( distance_2 < current_result.GetDistance2() ) {
		current_result.Set( *p_point );
		current_result.SetDistance2( distance_2 );
	      }
	      p_point++;
	    }
          } // for i_x
          y_position += mCells.GetNumberOfCells( 0 );
        } // for i_y
      } // for i_z
      radius *= 2.00;
    }
    // std::cout << " ... result = " << ( *current_result.Get()) << std::endl;
    // exit( 1);
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
      for ( std::size_t offset = mCells.GetCellBeginIndex( index );
            offset < mCells.GetCellEndIndex( index ); offset++ ) {
        if ( mpPoints[ offset ] == nullptr ) {
          mpPoints[ offset ] = &( *i_point );
          break;
        }
      }
    }
  }

  void SearchNearestInCell( std::size_t CellIndex, TObjectType const &ThePoint,
                            ResultType &rCurrentResult ) {
    for ( std::size_t offset = mCells.GetCellBeginIndex( CellIndex );
          offset < mCells.GetCellEndIndex( CellIndex ); offset++ ) {
      TObjectType *p_point = mpPoints[ offset ];
      auto distance_2 = Distance2( *p_point, ThePoint );
      if ( distance_2 <= rCurrentResult.GetDistance2() ) {
        rCurrentResult.Set( p_point );
        rCurrentResult.SetDistance2( distance_2 );
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
