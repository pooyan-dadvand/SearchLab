/* -*- c++ -*- */
#pragma once

#include <unordered_set>
#include "bounding_box.h"
#include "interval_count.h"

class BinsCellsContainer {
  static constexpr int Dimension = 3;

protected:
  std::array< std::size_t, Dimension > mNumberOfCells;
  using InternalPointType = std::array< double, Dimension >;
  BoundingBox< InternalPointType > mBoundingBox;
  InternalPointType mCellSize;
  InternalPointType mInverseOfCellSize;
  std::vector< std::size_t > mCellsBeginIndices;
  
public:

  template < typename TIteratorType >
  BinsCellsContainer( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd,
		      const std::size_t GridSize[ 3] ) // at the moment unused
    : mNumberOfCells( { { 1, 1, 1 } } ), mBoundingBox( PointsBegin, PointsEnd ) {
    
    std::size_t approximated_number_of_cells = std::distance( PointsBegin, PointsEnd );
    
    if ( approximated_number_of_cells == 0 )
      return;

    CalculateCellSize( approximated_number_of_cells, GridSize );

    InitializeCellsBeginIndices( PointsBegin, PointsEnd );
  }

  void SetNumberOfCells( std::array< std::size_t, Dimension > const &TheNumberOfCells ) {
    mNumberOfCells = TheNumberOfCells;
  }

  void SetNumberOfCells( std::size_t Axis, std::size_t TheNumberOfCells ) {
    mNumberOfCells[ Axis ] = TheNumberOfCells;
  }

  void SetCellBeginIndex( std::size_t Index, std::size_t Value ) {
    mCellsBeginIndices[ Index ] = Value;
  }

  std::size_t GetNumberOfCells( std::size_t Axis ) const { return mNumberOfCells[ Axis ]; }

  std::size_t GetTotalNumberOfCells() const { return mCellsBeginIndices.size(); }

  double GetCellSize( std::size_t Axis ) const { return mCellSize[ Axis ]; }

  std::size_t GetCellBeginIndex( std::size_t Index ) const { 
      return mCellsBeginIndices[ Index ];
  }

  template < typename TPointType >
  std::size_t CalculateCellIndex( TPointType const &ThePoint ) const {
    std::size_t result = 0;
    for ( int i_dim = Dimension - 1; i_dim > 0; i_dim-- ) {
      result += CalculatePosition( ThePoint[ i_dim ], i_dim );
      result *= mNumberOfCells[ i_dim - 1 ];
    }
    result += CalculatePosition( ThePoint[ 0 ], 0 );
    return result;
  }


  std::size_t CalculatePosition( double Coordinate, int ThisDimension ) const {
    auto distance = Coordinate - mBoundingBox.GetMinPoint()[ ThisDimension ];
    distance = ( distance < 0.00 ) ? 0.00 : distance;
    std::size_t position =
        static_cast< std::size_t >( distance * mInverseOfCellSize[ ThisDimension ] );
    std::size_t result= ( position > mNumberOfCells[ ThisDimension ] - 1 )
                            ? mNumberOfCells[ ThisDimension ] - 1
                            : position;
    return result;
  }

  void PrintStatistics() const;
  void PrintStatisticsStdVector() const;
  
private:
  void CalculateCellSize( std::size_t ApproximatedSize, const std::size_t GridSize[ 3] ) {
    std::size_t average_number_of_cells = static_cast< std::size_t >(
        std::pow( static_cast< double >( ApproximatedSize ), 1.00 / Dimension ) );
    std::array< double, 3 > lengths;
    double average_length = 0.00;
    for ( int i = 0; i < Dimension; i++ ) {
      lengths[ i ] = mBoundingBox.GetMaxPoint()[ i ] - mBoundingBox.GetMinPoint()[ i ];
      average_length += lengths[ i ];
    }
    average_length *= 1.00 / 3.00;

    if ( ( GridSize[ 0] == 0) || ( GridSize[ 1] == 0) || ( GridSize[ 2] == 0)) {
      if ( average_length < std::numeric_limits< double >::epsilon() ) {
	SetNumberOfCells( { { 1, 1, 1 } } );
	return;
      }
      
      for ( int i = 0; i < Dimension; i++ ) {
	SetNumberOfCells( i, static_cast< std::size_t >( lengths[ i] / average_length * ( double)average_number_of_cells) + 1);
	if ( mNumberOfCells[ i ] > 1 )
	  mCellSize[ i ] = lengths[ i ] / ( double)mNumberOfCells[ i ];
	else
	  mCellSize[ i ] = average_length;
	
	mInverseOfCellSize[ i ] = 1.00 / mCellSize[ i ];
      }
    } else { // User defined grid size
      for ( int i = 0; i < Dimension; i++ ) {
	SetNumberOfCells( i, GridSize[ i]);
	if ( mNumberOfCells[ i ] > 1 )
	  mCellSize[ i ] = lengths[ i ] / ( double)mNumberOfCells[ i ];
	else
	  mCellSize[ i ] = average_length;
	
	mInverseOfCellSize[ i ] = 1.00 / mCellSize[ i ];
      }
    }
  }

protected:
  template < typename TIteratorType >
  void InitializeCellsBeginIndices( TIteratorType const &PointsBegin,
                                    TIteratorType const &PointsEnd ) {
    mCellsBeginIndices.resize( mNumberOfCells[ 0 ] * mNumberOfCells[ 1 ] * mNumberOfCells[ 2 ] + 1, 0 );
    for ( auto i_point = PointsBegin; i_point != PointsEnd; i_point++ ) {
      mCellsBeginIndices[ CalculateCellIndex( *i_point ) + 1 ]++;
    }
    
    for ( std::size_t i_cell_begin = 1; i_cell_begin < mCellsBeginIndices.size(); i_cell_begin++ ) {
      mCellsBeginIndices[ i_cell_begin ] += mCellsBeginIndices[ i_cell_begin - 1 ];
    }
  }
};

inline void BinsCellsContainer::PrintStatisticsStdVector() const {
    std::size_t numUsedCells = 0;
    std::size_t lastOffset = 0;
    std::size_t totNumPoints = 0;
    std::size_t minNumPoints = this->GetTotalNumberOfCells();
    std::size_t maxNumPoints = 0;
    std::size_t numCellsWithSinglePoint = 0;
    for ( std::size_t idx = 1; idx < this->GetTotalNumberOfCells(); idx++ ) {
      lastOffset = this->GetCellBeginIndex( idx );
      std::size_t numberOfPoints = lastOffset - this->GetCellBeginIndex( idx - 1 );
      if ( numberOfPoints != 0 ) {
        numUsedCells++;
        totNumPoints += numberOfPoints;
        minNumPoints = ( numberOfPoints < minNumPoints ) ? numberOfPoints : minNumPoints;
        maxNumPoints = ( numberOfPoints > maxNumPoints ) ? numberOfPoints : maxNumPoints;
        if ( numberOfPoints == 1 )
          numCellsWithSinglePoint++;
      }
    }
    // the last this->GetTotalNumberOfCells() is already the total number of points...
    
    double occupancy_percent =
        ( 100.0 * ( ( double )numUsedCells / ( double )this->GetTotalNumberOfCells() ) );
    std::cout << " with " << numUsedCells << " used cells = " << occupancy_percent << " % occupancy"
              << std::endl;
    double bins_cells_array_size_MB =
        ( double )( this->GetTotalNumberOfCells() * sizeof( std::size_t ) ) / ( 1024.0 * 1024.0 );
    std::cout << " Bins size cell array = " << bins_cells_array_size_MB
              << " MB, used = " << bins_cells_array_size_MB * occupancy_percent / 100.0 << " MB"
              << std::endl;
    std::cout << " Number of points per cell ( Min, Avg, Max) = ( " << minNumPoints << ", "
              << ( double )totNumPoints / ( double )( numUsedCells ) << ", " << maxNumPoints << ")"
              << std::endl;

    std::cout << "Number of cells with only one point = " << numCellsWithSinglePoint << std::endl;
    IntervalCount ic( 8, ( double )minNumPoints, ( double )maxNumPoints );
    for ( std::size_t idx = 1; idx < this->GetTotalNumberOfCells(); idx++ ) {
      lastOffset = this->GetCellBeginIndex( idx );
      std::size_t numberOfPoints2 = lastOffset - this->GetCellBeginIndex( idx - 1 );
      if ( numberOfPoints2 != 0 ) {
        ic.countSample( ( double )numberOfPoints2 );
      }
    }
    ic.print();
}

inline void BinsCellsContainer::PrintStatistics() const {
    // Bins statistics
    std::cout << "=== Bins statistics === \n";
    std::locale prev_loc = std::cout.getloc();
    std::cout.imbue( std::locale( "" ) ); // for thousand separators ...
    std::cout << "Bin of ";
    std::size_t numberOfCells = 1;
    for ( std::size_t i = 0; i < Dimension; i++ ) {
      numberOfCells *= this->GetNumberOfCells( i );
      std::cout << this->GetNumberOfCells( i );
      if ( i < Dimension - 1 )
        std::cout << " x ";
    }
    std::cout << " = " << numberOfCells << " cells" << std::endl;
    std::cout << " = " << this->GetTotalNumberOfCells() << " cells" << std::endl;

    std::cout << "Using std::vector" << std::endl;
    this->PrintStatisticsStdVector();
    
    std::cout.imbue( prev_loc ); // restore previous locale, i.e. without thousand separators
    std::cout << "=== End of statistics === \n";
}
