/* -*- c++ -*- */
#pragma once

#include <iostream>

template < std::size_t dim_type > class SphereObject {

public:
  double coord[ dim_type ];
  double radius;
  std::size_t id;
  std::size_t tag;
  // int id;

  double &operator[]( std::size_t i ) { return coord[ i ]; }

  double const &operator[]( std::size_t i ) const { return coord[ i ]; }

  void operator=( SphereObject< dim_type > const &Other ) {
    for ( std::size_t i = 0; i < dim_type; i++ )
      coord[ i ] = Other.coord[ i ];
    radius = Other.radius;
  }

  SphereObject &Coordinates() { return *this; }

  SphereObject const &Coordinates() const { return *this; }
};

template < std::size_t dim_type >
std::ostream &operator<<( std::ostream &rOut, SphereObject< dim_type > &rObject ) {
  rOut << "(" << rObject.id << ") ";
  for ( std::size_t i = 0; i < dim_type; i++ )
    rOut << rObject[ i ] << " ";
  rOut << rObject.radius << " ";
  return rOut;
}

template < std::size_t dim_type >
std::istream &operator>>( std::istream &rIn, SphereObject< dim_type > &rObject ) {
  for ( std::size_t i = 0; i < dim_type; i++ )
    rIn >> rObject[ i ];
  rIn >> rObject.radius;

  return rIn;
}
