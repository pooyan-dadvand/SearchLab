/* -*- c++ -*- */
#pragma once

#include <iostream>

class Point {

public:
  double coord[ 3 ];
  std::size_t id;
  std::size_t tag;
  // int id;

  Point() {}
  Point( double x, double y, double z ) {
    coord[ 0 ] = x;
    coord[ 1 ] = y;
    coord[ 2 ] = z;
  }

  double &operator[]( std::size_t i ) { return coord[ i ]; }

  double const &operator[]( std::size_t i ) const { return coord[ i ]; }

  void operator=( Point const &Other ) {
    for ( std::size_t i = 0; i < 3; i++ )
      coord[ i ] = Other.coord[ i ];
  }

  Point operator+( const Point &t ) const {
    Point ret( *this );
    ret[ 0 ] += t[ 0 ];
    ret[ 1 ] += t[ 1 ];
    ret[ 2 ] += t[ 2 ];
    return ret;
  }

  Point operator-( const Point &t ) const {
    Point ret( *this );
    ret[ 0 ] -= t[ 0 ];
    ret[ 1 ] -= t[ 1 ];
    ret[ 2 ] -= t[ 2 ];
    return ret;
  }

  Point &Coordinates() { return *this; }

  Point const &Coordinates() const { return *this; }
};

std::ostream &operator<<( std::ostream &rOut, const Point &rPoint ) {
  rOut << "(" << rPoint.id << ") ";
  for ( std::size_t i = 0; i < 3; i++ )
    rOut << rPoint[ i ] << " ";
  return rOut;
}

std::istream &operator>>( std::istream &rIn, Point &rPoint ) {
  for ( std::size_t i = 0; i < 3; i++ )
    rIn >> rPoint[ i ];

  return rIn;
}
