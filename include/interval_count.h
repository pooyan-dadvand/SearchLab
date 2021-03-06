/* -*- c++ -*- */
#pragma once

#include <assert.h>
#include <iostream>
#include <fcntl.h>

class IntervalCount {
public:
  IntervalCount(): m_pivots( NULL), m_counts( NULL) {}
  IntervalCount( int numIntervals, double minPivot, double maxPivot): 
    m_pivots( NULL), m_counts( NULL), m_accumulated_sample_values( 0.0) {
    defineIntervals( numIntervals, minPivot, maxPivot);
  }
  IntervalCount( int numIntervals, const double *lstPivots): 
    m_pivots( NULL), m_counts( NULL), m_accumulated_sample_values( 0.0) {
    defineIntervals( numIntervals, lstPivots);
  }
  void defineIntervals( int numIntervals, double minPivot, double maxPivot);
  void defineIntervals( int numIntervals, const double *lstPivots);
  void countSample( double value);
  void print() const;
  void printAsFile( FILE *fo) const;
private:
  int m_numIntervals;
  double *m_pivots;
  int *m_counts;
  double m_accumulated_sample_values;
};

void IntervalCount::defineIntervals( int numIntervals, double minPivot, double maxPivot) {
  m_numIntervals = numIntervals;
  if ( m_pivots) { delete[] m_pivots; m_pivots = nullptr;}
  if ( m_counts) { delete[] m_counts; m_counts = nullptr;}
  if ( m_numIntervals) {
    double delta = ( maxPivot - minPivot) / ( double)m_numIntervals;
    if ( delta == 0.0) {
      m_numIntervals = 1; // all samples in the same interval ...
    }
    m_pivots = new double[ m_numIntervals + 1];
    m_counts = new int[ m_numIntervals];
    for ( int i = 0; i < m_numIntervals; i++) {
      m_pivots[ i] = ( double)( int)( 0.5 + ( double)i * delta + minPivot);
      m_counts[ i] = 0;
    }
    m_pivots[ m_numIntervals] = maxPivot;
  } // if ( m_numIntervals)
}

void IntervalCount::defineIntervals( int numIntervals, const double *lstPivots) {
  m_numIntervals = numIntervals;
  if ( m_pivots) { delete[] m_pivots; m_pivots = nullptr;}
  if ( m_counts) { delete[] m_counts; m_counts = nullptr;}
  if ( m_numIntervals) {
    m_pivots = new double[ m_numIntervals + 1];
    m_counts = new int[ m_numIntervals];
    for ( int i = 0; i < m_numIntervals; i++) {
      m_pivots[ i] = lstPivots[ i];
      m_counts[ i] = 0;
    }
    m_pivots[ m_numIntervals] = lstPivots[ m_numIntervals];
  } // if ( m_numIntervals)
}

void IntervalCount::countSample( double value) {
  assert( m_counts);
  assert( m_pivots);
  assert( ( value >= m_pivots[ 0]) && ( value <= m_pivots[ m_numIntervals]));
  if ( m_numIntervals) {
    for ( int i = 1; i < m_numIntervals + 1; i++) {
      if ( value <= m_pivots[ i]) {
	m_counts[ i - 1]++; // 1 cell with that many points
	m_accumulated_sample_values += value; // add number of points
	break;
      }
    }
  }
}

void IntervalCount::print() const {
  std::cout << "Num intervals = " << m_numIntervals << std::endl;
  if ( m_numIntervals) {
    std::cout << "Pivots";
    for ( int i = 0; i < m_numIntervals + 1; i++) {
      std::cout << "\t" << ( int)m_pivots[ i];
    }
    std::cout << std::endl;
    std::cout << "Counts";
    int totalCounts = 0;
    for ( int i = 0; i < m_numIntervals; i++) {
      std::cout << "\t" << m_counts[ i];
      totalCounts += m_counts[ i];
    }
    std::cout << std::endl;
    std::cout << "Total counts = " << totalCounts << std::endl;
  }
}

void IntervalCount::printAsFile( FILE *fo) const {
  if ( m_numIntervals) {
    std::size_t total_counts = 0;
    for ( int i = 0; i < m_numIntervals; i++) {
      total_counts +=m_counts[ i];
    }
    double factor_points = ( double)GetTotalNumberOfPoints() / m_accumulated_sample_values;
    double factor_cells = ( double)GetMaxNumberOfCells()/ ( double)GetCurrentNumberOfCells();
    fprintf( fo, "# Number of points = %lld\n", ( long long int)m_accumulated_sample_values);
    fprintf( fo, "#_points/cell   count   #_points/cell_normalized   count_normalized\n");
    for ( int i = 0; i < m_numIntervals; i++) {
      fprintf( fo, "%d   %d   %g   %g\n", ( int)m_pivots[ i + 1], m_counts[ i],
	       ( double)( int)m_pivots[ i + 1] * factor_points, // 1st to int to get the integer part only
	       ( double)m_counts[ i] * factor_cells);
    }
  }
  fprintf( fo, "\n\n");
}
