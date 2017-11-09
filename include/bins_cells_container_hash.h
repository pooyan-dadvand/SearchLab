/* -*- c++ -*- */
#pragma once

#include <unordered_set>
#include "bounding_box.h"
#include "interval_count.h"
#include "parallel_coherent_hash.h"

class BinsCellsContainerHash {
  static constexpr int Dimension = 3;

protected:
  std::array< std::size_t, Dimension > mNumberOfCells;
  using InternalPointType = std::array< double, Dimension >;
  BoundingBox< InternalPointType > mBoundingBox;
  InternalPointType mCellSize;
  InternalPointType mInverseOfCellSize;
  ParallelCoherentHash< std::size_t, std::size_t> m_PCHCellsBeginIndices;
  std::size_t m_numCells;
  
public:
  template < typename TIteratorType >
  BinsCellsContainerHash( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd )
    : mNumberOfCells( { { 1, 1, 1 } } ), mBoundingBox( PointsBegin, PointsEnd ), m_numCells( 0) {
    
    std::size_t approximated_number_of_cells = std::distance( PointsBegin, PointsEnd );
    
    if ( approximated_number_of_cells == 0 )
      return;

    CalculateCellSize( approximated_number_of_cells );

    InitializeCellsBeginIndices( PointsBegin, PointsEnd );
  }

  void SetNumberOfCells( std::array< std::size_t, Dimension > const &TheNumberOfCells ) {
    mNumberOfCells = TheNumberOfCells;
  }

  void SetNumberOfCells( std::size_t Axis, std::size_t TheNumberOfCells ) {
    mNumberOfCells[ Axis ] = TheNumberOfCells;
  }

  void SetCellBeginIndex( std::size_t Index, std::size_t Value ) {
    m_PCHCellsBeginIndices.getDataRef( Index ) = Value;
  }

  std::size_t GetNumberOfCells( std::size_t Axis ) const { return mNumberOfCells[ Axis ]; }

  std::size_t GetTotalNumberOfCells() const { return m_numCells;}//mCellsBeginIndices.size(); }

  double GetCellSize( std::size_t Axis ) const { return mCellSize[ Axis ]; }

  std::size_t GetCellBeginIndex( std::size_t Index ) const { 
      return m_PCHCellsBeginIndices.getData( Index);
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
    return ( position > mNumberOfCells[ ThisDimension ] - 1 ) ? mNumberOfCells[ ThisDimension ] - 1
                                                              : position;
  }

  void PrintStatistics() const;
  void PrintStatisticsHash() const;
  
private:
  void CalculateCellSize( std::size_t ApproximatedSize ) {
    std::size_t avarage_number_of_cells = static_cast< std::size_t >(
        std::pow( static_cast< double >( ApproximatedSize ), 1.00 / Dimension ) );
    std::array< double, 3 > lengths;
    double avarage_length = 0.00;
    for ( int i = 0; i < Dimension; i++ ) {
      lengths[ i ] = mBoundingBox.GetMaxPoint()[ i ] - mBoundingBox.GetMinPoint()[ i ];
      avarage_length += lengths[ i ];
    }
    avarage_length *= 1.00 / 3.00;

    if ( avarage_length < std::numeric_limits< double >::epsilon() ) {
      SetNumberOfCells( { { 1, 1, 1 } } );
      return;
    }

    for ( int i = 0; i < Dimension; i++ ) {
      SetNumberOfCells( i, static_cast< std::size_t >( lengths[ i] / avarage_length * ( double)avarage_number_of_cells) + 1);
      if ( mNumberOfCells[ i ] > 1 )
        mCellSize[ i ] = lengths[ i ] / ( double)mNumberOfCells[ i ];
      else
        mCellSize[ i ] = avarage_length;

      mInverseOfCellSize[ i ] = 1.00 / mCellSize[ i ];
    }
  }

protected:
  template < typename TIteratorType >
  void InitializeCellsBeginIndices( TIteratorType const &PointsBegin,
                                    TIteratorType const &PointsEnd ) {
    m_numCells = mNumberOfCells[ 0 ] * mNumberOfCells[ 1 ] * mNumberOfCells[ 2 ] + 1;
    
    // first we need to calculate the number of used cells
    std::unordered_set< size_t> setUsedIndices;
    for ( auto i_point = PointsBegin; i_point != PointsEnd; i_point++ ) {
      std::size_t idx = CalculateCellIndex( *i_point ) + 1;
      setUsedIndices.insert( idx);
    }
    m_numCells = setUsedIndices.size(); // number of unique cells
    m_PCHCellsBeginIndices.resize( m_numCells, 0);
    for ( auto i_point = PointsBegin; i_point != PointsEnd; i_point++ ) {
      std::size_t idx = CalculateCellIndex( *i_point ) + 1;
      m_PCHCellsBeginIndices.getDataRef( idx )++;
    }
    for ( std::size_t idx = 1; idx < m_PCHCellsBeginIndices.getRawHashTableSize(); idx++ ) {
      m_PCHCellsBeginIndices.getRawEntryDataRef(  idx) += m_PCHCellsBeginIndices.getRawEntryData( idx - 1);
    }
  }
};

inline void BinsCellsContainerHash::PrintStatisticsHash() const {
    m_PCHCellsBeginIndices.PrintStatistics();
    std::size_t numUsedCells = 0;
    std::size_t lastOffset = 0;
    std::size_t totNumPoints = 0;
    std::size_t minNumPoints = this->GetTotalNumberOfCells(); // a big number like any other...
    std::size_t maxNumPoints = 0;
    std::size_t numCellsWithSinglePoint = 0;
    std::size_t last_idx = 0;
    for ( std::size_t idx = 1; idx < m_PCHCellsBeginIndices.getRawHashTableSize(); idx++ ) {
      if ( m_PCHCellsBeginIndices.getRawEntryUsed( idx)) {
	lastOffset = m_PCHCellsBeginIndices.getRawEntryData( idx );
	std::cout << " entry " << idx << " = " << lastOffset << std::endl;
	std::size_t numberOfPoints = lastOffset - m_PCHCellsBeginIndices.getRawEntryData( last_idx );
	last_idx = idx;
	if ( numberOfPoints != 0 ) {
	  numUsedCells++;
	  totNumPoints += numberOfPoints;
	  minNumPoints = ( numberOfPoints < minNumPoints ) ? numberOfPoints : minNumPoints;
	  maxNumPoints = ( numberOfPoints > maxNumPoints ) ? numberOfPoints : maxNumPoints;
	  if ( numberOfPoints == 1 )
	    numCellsWithSinglePoint++;
	}
      }
    }
    // the last this->GetTotalNumberOfCells() is already the total number of points...
    
    double occupancy_percent =
        ( 100.0 * ( ( double )numUsedCells / ( double )m_PCHCellsBeginIndices.getRawHashTableSize() ) );
    std::cout << " with " << numUsedCells << " used cells = " << occupancy_percent << " % occupancy"
              << std::endl;
    std::cout << " Number of points per cell ( Min, Avg, Max) = ( " << minNumPoints << ", "
              << ( double )totNumPoints / ( double )( numUsedCells ) << ", " << maxNumPoints << ")"
              << std::endl;

    std::cout << "Number of cells with only one point = " << numCellsWithSinglePoint << std::endl;
    IntervalCount ic( 8, ( double )minNumPoints, ( double )maxNumPoints );
    last_idx = 0;
    for ( std::size_t idx = 1; idx < m_PCHCellsBeginIndices.getRawHashTableSize(); idx++ ) {
      if ( m_PCHCellsBeginIndices.getRawEntryUsed( idx)) {
	lastOffset = m_PCHCellsBeginIndices.getRawEntryData( idx );
	std::size_t numberOfPoints = lastOffset - m_PCHCellsBeginIndices.getRawEntryData( last_idx );
	last_idx = idx;
	if ( numberOfPoints != 0 ) {
	  ic.countSample( ( double )numberOfPoints );
	}
      }
    }
    ic.print();
}

inline void BinsCellsContainerHash::PrintStatistics() const {
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

    // for ( std::size_t idx = 0; idx < Dimension; idx++) {
    //   std::cout << " Cell size in axis = " << idx << " is " << this->GetCellSize( idx) <<
    //   std::endl;
    // }

    std::cout << "Using ParallelCoherentHash" << std::endl;
    this->PrintStatisticsHash();

    std::cout.imbue( prev_loc ); // restore previous locale, i.e. without thousand separators
    std::cout << "=== End of statistics === \n";
  }
