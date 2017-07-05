#pragma once

#include <iostream>

template< std::size_t dim_type>
class Point {

public:

	double       coord[dim_type];
	std::size_t  id;
	std::size_t  tag;
	//int id;

  Point() {}
  Point(double x, double y, double z) {
    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
  }

	double& operator[](std::size_t i) { return coord[i]; }

	double const & operator[](std::size_t i) const { return coord[i]; }

	void operator=(Point<dim_type> const& Other) {
		for (std::size_t i = 0; i < dim_type; i++)
			coord[i] = Other.coord[i];
	}

  Point<dim_type> operator+ (const Point<dim_type>& t) const
  {
    Point<dim_type> ret(*this);
    ret[0] += t[0];
    ret[1] += t[1];
    ret[2] += t[2];
    return ret;
  }

  Point<dim_type> operator- (const Point<dim_type>& t) const
  {
    Point<dim_type> ret(*this);
    ret[0] -= t[0];
    ret[1] -= t[1];
    ret[2] -= t[2];
    return ret;
  }

	Point& Coordinates() { return *this; }

	Point const& Coordinates() const { return *this; }
};

template< std::size_t dim_type >
std::ostream & operator<<(std::ostream& rOut, Point<dim_type> & rPoint) {
	rOut << "(" << rPoint.id << ") ";
	for (std::size_t i = 0; i < dim_type; i++)
		rOut << rPoint[i] << " ";
	return rOut;
};

template< std::size_t dim_type >
std::istream & operator>>(std::istream& rIn, Point<dim_type> & rPoint) {
	for (std::size_t i = 0; i < dim_type; i++)
		rIn >> rPoint[i];

	return rIn;
};
