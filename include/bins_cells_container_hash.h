/* -*- c++ -*- */
#pragma once

#include <unordered_set>
#include <algorithm>
#include "bounding_box.h"
#include "interval_count.h"
#include "parallel_coherent_hash.h"

class BinsCellsContainerHash {
  static constexpr int Dimension = 3;

  using t_Offset = std::size_t;
  typedef struct {
    t_Offset offset_ini, offset_end;
  } t_CellContents;
  
protected:
  std::array< std::size_t, Dimension > mNumberOfCells;
  using InternalPointType = std::array< double, Dimension >;
  BoundingBox< InternalPointType > mBoundingBox;
  InternalPointType mCellSize;
  InternalPointType mInverseOfCellSize;
  ParallelCoherentHash< t_CellContents, std::size_t> m_PCHCellsBeginIndices;
  std::size_t m_numCells; // = # cells of the full grid
  std::size_t m_numUsedCells;
  std::vector< std::size_t> m_lstUsedIndices;

public:
  template < typename TIteratorType >
  BinsCellsContainerHash( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd )
    : mNumberOfCells( { { 1, 1, 1 } } ), mBoundingBox( PointsBegin, PointsEnd ), 
      m_numCells( 0), m_numUsedCells( 0) {
    CalculateBinsCells( PointsBegin, PointsEnd);
  }

  template < typename TIteratorType >
  BinsCellsContainerHash( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd,
			  const BoundingBox< InternalPointType > &BBox )
    : mNumberOfCells( { { 1, 1, 1 } } ), mBoundingBox( BBox ), 
      m_numCells( 0), m_numUsedCells( 0) {
    CalculateBinsCells( PointsBegin, PointsEnd);
  }

  // to specify a custom gridSize
  template < typename TIteratorType >
  BinsCellsContainerHash( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd,
			  const std::size_t GridSize[ 3] )
    : mNumberOfCells( { { 1, 1, 1 } } ), mBoundingBox( PointsBegin, PointsEnd ), 
      m_numCells( 0), m_numUsedCells( 0) {
    CalculateBinsCellsForUserGridSize( PointsBegin, PointsEnd, GridSize);
  }
  template < typename TIteratorType >
  BinsCellsContainerHash( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd,
			  const std::size_t GridSize[ 3], 
			  const BoundingBox< InternalPointType > &BBox )
    : mNumberOfCells( { { 1, 1, 1 } } ), mBoundingBox( BBox ), 
      m_numCells( 0), m_numUsedCells( 0) {
    CalculateBinsCellsForUserGridSize( PointsBegin, PointsEnd, GridSize);
  }

  std::size_t GetNumberOfCells( std::size_t Axis ) const { return mNumberOfCells[ Axis ]; }
  std::size_t GetTotalNumberOfCells() const { return m_numCells;}
  std::size_t GetNumberOfUsedCells() const { return m_numUsedCells;}
  double GetCellSize( std::size_t Axis ) const { return mCellSize[ Axis ]; }
  BoundingBox< InternalPointType > GetBoundingBox() const { return mBoundingBox;}

  // std::size_t GetCellBeginIndex( std::size_t Index ) const {
  //   bool found = false;
  //   std::size_t ret = m_PCHCellsBeginIndices.getData( Index, found).offset_ini;
  //   if ( !found)
  //     ret = ( std::size_t)-1;
  //   return ret;
  // }
  // 
  // std::size_t GetCellEndIndex( std::size_t Index ) const { 
  //   bool found = false;
  //   std::size_t ret = m_PCHCellsBeginIndices.getData( Index, found).offset_end;
  //   if ( !found)
  //     ret = ( std::size_t)-1;
  //   return ret;
  // }

  bool GetCellStoredOffsets( std::size_t Index, std::size_t &begin, std::size_t &end) const {
    bool found = false;
    t_CellContents ret = m_PCHCellsBeginIndices.getData( Index, found);
    if ( !found) {
      return false;
    }
    begin = ret.offset_ini;
    end = ret.offset_end;
    return true;
  }

  const std::vector< std::size_t> &GetListUsedCellIndices() const { return m_lstUsedIndices;}

  bool CellIsEmpty( std::size_t Index ) const {
    bool found = false;
    m_PCHCellsBeginIndices.getData( Index, found);
    return found;
  }
  bool CellIsNotEmpty( std::size_t Index ) const { return !CellIsEmpty( Index);}

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

  InternalPointType CalculateCentreOfCell( std::size_t Index ) const {
    const std::size_t num_xy = mNumberOfCells[ 0] * mNumberOfCells[ 1];
    const std::size_t idx_z = Index / num_xy;
    const std::size_t rem =  Index % num_xy;
    const std::size_t idx_y = rem / mNumberOfCells[ 0];
    const std::size_t idx_x = rem % mNumberOfCells[ 0];
    const InternalPointType &min = mBoundingBox.GetMinPoint();
    InternalPointType centre = { { 
	( double)idx_x * mCellSize[ 0] + min[ 0], 
	( double)idx_y * mCellSize[ 1] + min[ 1], 
	( double)idx_z * mCellSize[ 2] + min[ 2]}};
    return centre;
  }

  void GetCellVectorIndices( std::size_t Index, std::size_t &idx_x, std::size_t &idx_y, std::size_t &idx_z ) {
    const std::size_t num_xy = mNumberOfCells[ 0 ] * mNumberOfCells[ 1 ];
    idx_z = Index / num_xy;
    const std::size_t rem = Index % num_xy;
    idx_y = rem / mNumberOfCells[ 0 ];
    idx_x = rem % mNumberOfCells[ 0 ];
  }

  std::size_t CalculatePosition( double Coordinate, int ThisDimension ) const {
    auto distance = Coordinate - mBoundingBox.GetMinPoint()[ ThisDimension ];
    distance = ( distance < 0.00 ) ? 0.00 : distance;
    std::size_t position =
        static_cast< std::size_t >( distance * mInverseOfCellSize[ ThisDimension ] );
    std::size_t result = ( position > mNumberOfCells[ ThisDimension ] - 1 )
                             ? mNumberOfCells[ ThisDimension ] - 1
                             : position;
    return result;
  }

  void PrintStatistics() const;
  void PrintStatisticsHash() const;
  void PrintDensitiesInFile( const char *filename) const;
  
private:

  template < typename TIteratorType >
  void CalculateBinsCells( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd) {
    std::size_t approximated_number_of_cells = std::distance( PointsBegin, PointsEnd );
    
    if ( approximated_number_of_cells == 0 )
      return;
    
    CalculateCellSize( approximated_number_of_cells );
    
    InitializeCellsBeginIndices( PointsBegin, PointsEnd );
  }

  template < typename TIteratorType >
  void CalculateBinsCellsForUserGridSize( TIteratorType const &PointsBegin, TIteratorType const &PointsEnd,
					  const std::size_t GridSize[ 3]) {
    std::size_t approximated_number_of_cells = std::distance( PointsBegin, PointsEnd );
    
    if ( approximated_number_of_cells == 0 )
      return;
    
    CalculateCellSizeForUserGridSize( GridSize );
    
    InitializeCellsBeginIndices( PointsBegin, PointsEnd );
  }

  void SetNumberOfCells( std::array< std::size_t, Dimension > const &TheNumberOfCells ) {
    mNumberOfCells = TheNumberOfCells;
  }
  void SetNumberOfCells( std::size_t Axis, std::size_t TheNumberOfCells ) {
    mNumberOfCells[ Axis ] = TheNumberOfCells;
  }
  // void SetCellBeginIndex( std::size_t Index, std::size_t Value ) {
  //   m_PCHCellsBeginIndices.getDataRef( Index ) = Value;
  // }
  void CalculateCellSize( std::size_t ApproximatedSize) {
    double average_length = 0.00;
    std::array< double, 3 > lengths;
    for ( int i = 0; i < Dimension; i++ ) {
      lengths[ i ] = mBoundingBox.GetMaxPoint()[ i ] - mBoundingBox.GetMinPoint()[ i ];
      average_length += lengths[ i ];
    }
    average_length *= 1.00 / 3.00;

    std::size_t average_number_of_cells = 
      static_cast< std::size_t >( std::pow( static_cast< double >( ApproximatedSize ), 1.00 / Dimension ) );
    
    if ( average_length < std::numeric_limits< double >::epsilon() ) {
      SetNumberOfCells( { { 1, 1, 1 } } );
      return;
    }
    
    for ( int i = 0; i < Dimension; i++ ) {
      SetNumberOfCells( i, static_cast< std::size_t >( lengths[ i] / average_length * ( double)average_number_of_cells + 1));
      if ( mNumberOfCells[ i ] > 1 )
	mCellSize[ i ] = lengths[ i ] / ( double)mNumberOfCells[ i ];
      else
	mCellSize[ i ] = average_length;
      
      mInverseOfCellSize[ i ] = 1.00 / mCellSize[ i ];
    }  
  }

  void CalculateCellSizeForUserGridSize( const std::size_t GridSize[ 3] ) {
    double average_length = 0.00;
    std::array< double, 3 > lengths;
    for ( int i = 0; i < Dimension; i++ ) {
      lengths[ i ] = mBoundingBox.GetMaxPoint()[ i ] - mBoundingBox.GetMinPoint()[ i ];
      average_length += lengths[ i ];
    }
    average_length *= 1.00 / 3.00;
    
    // User defined grid size
    for ( int i = 0; i < Dimension; i++ ) {
      SetNumberOfCells( i, GridSize[ i]);
      if ( mNumberOfCells[ i ] > 1 )
	mCellSize[ i ] = lengths[ i ] / ( double)mNumberOfCells[ i ];
      else
	mCellSize[ i ] = average_length;
      
      mInverseOfCellSize[ i ] = 1.00 / mCellSize[ i ];
    }
  }

protected:
  template < typename TIteratorType >
  void InitializeCellsBeginIndices( TIteratorType const &PointsBegin,
                                    TIteratorType const &PointsEnd );
};


template < typename TIteratorType >
inline void BinsCellsContainerHash::InitializeCellsBeginIndices( TIteratorType const &PointsBegin,
								 TIteratorType const &PointsEnd ) {
  m_numCells = mNumberOfCells[ 0 ] * mNumberOfCells[ 1 ] * mNumberOfCells[ 2 ];
  
  // first we need to calculate the number of used cells
  std::unordered_set< size_t> setUsedIndices;
  for ( auto i_point = PointsBegin; i_point != PointsEnd; i_point++ ) {
    std::size_t idx = CalculateCellIndex( *i_point );
    setUsedIndices.insert( idx);
  }
  m_numUsedCells = setUsedIndices.size(); // number of unique cells
  m_PCHCellsBeginIndices.resize( m_numUsedCells, { 0, 0});
  
  for ( auto i_point = PointsBegin; i_point != PointsEnd; i_point++ ) {
    std::size_t idx = CalculateCellIndex( *i_point );
    m_PCHCellsBeginIndices.getDataRef( idx ).offset_end++;
  }
  
  // sort used indices, to avoid loop over all cells: too costly for big grids
  m_lstUsedIndices.clear();
  for ( auto it_idx = setUsedIndices.begin(); it_idx != setUsedIndices.end(); it_idx++) {
    m_lstUsedIndices.push_back( *it_idx);
  }
  std::sort( m_lstUsedIndices.begin(), m_lstUsedIndices.end());
  
  bool first_time = true;
  std::size_t last_idx = 0;
  // for ( std::size_t idx = 0; idx < m_numCells; idx++) {
  for ( auto it_idx = m_lstUsedIndices.begin(); it_idx < m_lstUsedIndices.end(); it_idx++) {
    std::size_t idx = *it_idx;
    bool found = false;
    m_PCHCellsBeginIndices.getData( idx, found); // look if cell is there
    if ( !found) continue; // cell not stored 
    if ( first_time) {
      last_idx = idx;
      first_time = false;
      continue;
    }
    const std::size_t last_offset_end = m_PCHCellsBeginIndices.getData( last_idx).offset_end;
    m_PCHCellsBeginIndices.getDataRef( idx).offset_ini = last_offset_end;
    m_PCHCellsBeginIndices.getDataRef( idx).offset_end += last_offset_end;
    last_idx = idx;
  }
}

inline void BinsCellsContainerHash::PrintStatisticsHash() const {
    m_PCHCellsBeginIndices.PrintStatistics();
    std::size_t numUsedCells = 0;
    std::size_t totNumPoints = 0;
    std::size_t minNumPoints = this->GetTotalNumberOfCells(); // a big number like any other...
    std::size_t maxNumPoints = 0;
    std::size_t numCellsWithSinglePoint = 0;
    std::size_t numCellsWithLessThan10 = 0;
    std::size_t numCellsWithLessThan100 = 0;
    for ( std::size_t idx = 0; idx < m_PCHCellsBeginIndices.getRawHashTableSize(); idx++ ) {
      if ( m_PCHCellsBeginIndices.getRawEntryUsed( idx)) {
	std::size_t lastOffset = m_PCHCellsBeginIndices.getRawEntryData( idx ).offset_ini;
	std::size_t numberOfPoints = m_PCHCellsBeginIndices.getRawEntryData( idx ).offset_end - lastOffset;
	if ( numberOfPoints != 0 ) {
	  numUsedCells++;
	  totNumPoints += numberOfPoints;
	  minNumPoints = ( numberOfPoints < minNumPoints ) ? numberOfPoints : minNumPoints;
	  maxNumPoints = ( numberOfPoints > maxNumPoints ) ? numberOfPoints : maxNumPoints;
	  if ( numberOfPoints == 1 )
	    numCellsWithSinglePoint++;
          if ( numberOfPoints <= 10 )
            numCellsWithLessThan10++;
          if ( numberOfPoints <= 100 )
            numCellsWithLessThan100++;
	}
      }
    }
    
    double occupancy_percent =
        ( 100.0 * ( ( double )numUsedCells / ( double )m_PCHCellsBeginIndices.getRawHashTableSize() ) );
    std::cout << " with " << numUsedCells << " used cells = " << occupancy_percent << " % occupancy"
              << std::endl;
    std::cout << " Number of points per cell ( Min, Avg, Max) = ( " << minNumPoints << ", "
              << ( double )totNumPoints / ( double )( numUsedCells ) << ", " << maxNumPoints << ")"
              << std::endl;

    bool detailed_statistics = false; // true;
    if ( detailed_statistics) {
      std::cout << "Number of cells with only ( 1, <= 10, <= 100) points = ( "
                << numCellsWithSinglePoint << ", " 
                << numCellsWithLessThan10 << ", " 
                << numCellsWithLessThan100 << ") " << std::endl;
      
      IntervalCount ic( 8, ( double )minNumPoints, ( double )maxNumPoints );
      for ( std::size_t idx = 0; idx < m_PCHCellsBeginIndices.getRawHashTableSize(); idx++ ) {
        if ( m_PCHCellsBeginIndices.getRawEntryUsed( idx)) {
          std::size_t lastOffset = m_PCHCellsBeginIndices.getRawEntryData( idx ).offset_ini;
          std::size_t numberOfPoints = m_PCHCellsBeginIndices.getRawEntryData( idx ).offset_end - lastOffset;
          if ( numberOfPoints != 0 ) {
            ic.countSample( ( double )numberOfPoints );
          }
        }
      }
      ic.print();
    }
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

    std::cout << "Using ParallelCoherentHash" << std::endl;
    this->PrintStatisticsHash();

    std::cout.imbue( prev_loc ); // restore previous locale, i.e. without thousand separators
    std::cout << "=== End of statistics === \n";
  }

void BinsCellsContainerHash::PrintDensitiesInFile( const char *filename) const {
  double lstPivots[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
			 10, 20, 30, 40, 50, 60, 70, 80, 90, 
			 100, 200, 300, 400, 500, 600, 700, 800, 900, 
			 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 
			 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 
			 100000};

  if ( filename && *filename) {
    FILE *fo = fopen( filename, "a");
    if ( fo) {  
      fprintf( fo, "# Bin of %d x %d x %d = %lld cells.\n", 
	       ( int)( this->GetNumberOfCells( 0 )), ( int)( this->GetNumberOfCells( 1 )), ( int)( this->GetNumberOfCells( 2 )), 
	       ( long long)( this->GetNumberOfCells( 0 ) * this->GetNumberOfCells( 1 ) * this->GetNumberOfCells( 2 )));
      IntervalCount ic( ( int)( sizeof( lstPivots) / sizeof( double)) - 1, lstPivots); // num Intervals = num Pivots - 1
      for ( std::size_t idx = 0; idx < m_PCHCellsBeginIndices.getRawHashTableSize(); idx++ ) {
        if ( m_PCHCellsBeginIndices.getRawEntryUsed( idx)) {
          std::size_t lastOffset = m_PCHCellsBeginIndices.getRawEntryData( idx ).offset_ini;
          std::size_t numberOfPoints = m_PCHCellsBeginIndices.getRawEntryData( idx ).offset_end - lastOffset;
          if ( numberOfPoints != 0 ) {
            ic.countSample( ( double )numberOfPoints );
          }
        }
      }
      ic.printAsFile( fo);
      fclose( fo);
    }
  } else {
    std::cout << "Nothing to do in PrintDensitiesInFile, as filename is nullptr." << std::endl;
  }
}
