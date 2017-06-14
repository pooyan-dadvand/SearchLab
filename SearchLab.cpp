
//System includes
#include <fstream>
#include <iostream>
#include <ctime>
#include <string>
#include <cstdlib>

//Kratos Independent
#define KRATOS_INDEPENDENT

//Kratos includes
#include "spatial_containers/spatial_containers.h"

#include "points_bins.h"

double GetCurrentTime()
{
#ifndef _OPENMP
	return std::clock() / static_cast<double>(CLOCKS_PER_SEC);
#else
	return  omp_get_wtime();
#endif
}


double rrandom() {
	return double(rand()) / RAND_MAX;
};

template< std::size_t dim_type>
class Point {

public:

	double       coord[dim_type];
	std::size_t  id;
	std::size_t  tag;
	//int id;

	double& operator[](std::size_t i) { return coord[i]; }

	double const & operator[](std::size_t i) const { return coord[i]; }

	void RandomCoord() {
		for (std::size_t i = 0; i < dim_type; i++)
			coord[i] = rrandom();
	}

	void operator=(Point<dim_type> const& Other) {
		for (std::size_t i = 0; i < dim_type; i++)
			coord[i] = Other.coord[i];
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

template< class T, std::size_t dim >
class PointDistance {
public:
	double operator()(T const& p1, T const& p2) {
		double dist = 0.0;
		for (std::size_t i = 0; i < dim; i++) {
			double tmp = p1[i] - p2[i];
			dist += tmp*tmp;
		}
		return sqrt(dist);
	}
};

template< class T, std::size_t dim >
class PointDistance2 {
public:
	double operator()(T const& p1, T const& p2) {
		double dist = 0.0;
		for (std::size_t i = 0; i < dim; i++) {
			double tmp = p1[i] - p2[i];
			dist += tmp*tmp;
		}
		return dist;
	}
};

template< std::size_t dim >
bool LowerPoint(Point<dim> const& reference, Point<dim> const& new_) {
	for (std::size_t i = 0; i < dim; i++)
		if (reference[i] < new_[i])
			return false;
	return true;
};

template< std::size_t dim >
bool UpperPoint(Point<dim> const& reference, Point<dim> const& new_) {
	for (std::size_t i = 0; i < dim; i++)
		if (reference[i] > new_[i])
			return false;
	return true;
};


template< class TreeType, class PointType, class IteratorType, class DistanceIterator>
void RunTestOMP(char const Title[], IteratorType PBegin, IteratorType PEnd, IteratorType results0, DistanceIterator distances0, std::size_t MaxResults, PointType* allPoints, double radius0, std::size_t numsearch, std::size_t bucket_size)
{
	double t0, t1;
	std::size_t n;

	std::size_t npoints = PEnd - PBegin;
	std::vector<std::size_t> numresArray(numsearch);
	PointType& point0 = allPoints[0];

	std::size_t numsearch_nearest = 100 * numsearch;

	t0 = GetCurrentTime();
	TreeType   nodes_tree(PBegin, PEnd, bucket_size);
	t1 = GetCurrentTime();
	std::cout << Title << "\t" << t1 - t0 << "\t";

	t0 = GetCurrentTime();
#pragma omp parallel for
	for (std::size_t i = 0; i < numsearch; i++)
		n = nodes_tree.SearchInRadius(point0, radius0, results0, distances0, MaxResults);
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < numsearch; i++)
		n = nodes_tree.SearchInRadius(point0, radius0, results0, distances0, MaxResults);
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	PointType** PNearestArray = new PointType*[numsearch];
	std::vector<double> distancesArrayNearest(numsearch);
	PointType* PNearest = PNearestArray[0];

	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < 100; i++)
		nodes_tree.SearchNearestPoint(allPoints, numsearch, PNearestArray, distancesArrayNearest);
	t1 = GetCurrentTime();

	std::cout << t1 - t0 << "\t";
	double distance = 0;


	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < numsearch_nearest; i++)
		PNearest = nodes_tree.SearchNearestPoint(point0, distance);
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	std::cout << n << "\t" << *PNearestArray[0] << std::endl;
};

int main(int arg, char* argv[])
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//  Defines and initialization                                                                      //
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	static const std::size_t Dim = 3;

	typedef Point<Dim> PointType;

	typedef PointType*                              PtrPointType;
	typedef PtrPointType*                           PointVector;
	typedef PtrPointType*                           PointIterator;

	typedef double*                                 DistanceVector;
	typedef double*                                 DistanceIterator;

	//typedef DistanceSpatialContainersConfigure      ConfigType;
	//typedef Kratos::OctreeBinaryCell<ConfigType>    CellType;

	typedef Kratos::Bucket<         Dim,
		PointType,
		PointVector,
		PtrPointType,
		PointIterator,
		DistanceIterator,
		PointDistance2<PointType, Dim>
	>                                               bucket_type;                //Bucket;

	typedef Kratos::Bins<           Dim,
		PointType,
		PointVector,
		PtrPointType,
		PointIterator,
		DistanceIterator,
		PointDistance2<PointType, Dim> >                 static_bins_type;           //StaticBins;

	typedef Kratos::BinsDynamic<    Dim,
		PointType,
		PointVector,
		PtrPointType,
		PointIterator,
		DistanceIterator,
		PointDistance2<PointType, Dim> >                 dynamic_bins_type;          //DynamicBins;

																					 //typedef Kratos::OctreeBinary<CellType>                                          octree_binary_type;         //OctreeBinary;

	typedef Kratos::Tree< Kratos::KDTreePartition<bucket_type> >                    kdtree_type;                //Kdtree;
	typedef Kratos::Tree< Kratos::KDTreePartitionAverageSplit<bucket_type> >        kdtree_average_split_type;  //Kdtree;
	typedef Kratos::Tree< Kratos::KDTreePartitionMidPointSplit<bucket_type> >       kdtree_midpoint_split_type; //Kdtree;
	typedef Kratos::Tree< Kratos::OCTreePartition<bucket_type> >                    octree_type;                //Octree;
	typedef Kratos::Tree< Kratos::KDTreePartitionMidPointSplit<static_bins_type> >  kdtree_static_bins_type;    //KdtreeBins;
	typedef Kratos::Tree< Kratos::OCTreePartition<static_bins_type> >               octree_static_bins_type;    //OctreeBins;
	typedef Kratos::Tree< Kratos::KDTreePartitionMidPointSplit<dynamic_bins_type> > kdtree_dynamic_bins_type;   //KdtreeBins;
	typedef Kratos::Tree< Kratos::OCTreePartition<dynamic_bins_type> >              octree_bins_type;           //OctreeBins;

																												//////////////////////////////////////////////////////////////////////////////////////////////////////
																												// Input data                                                                                       //
																												//////////////////////////////////////////////////////////////////////////////////////////////////////

																												// set format
	std::cout.precision(4);

	PointVector points;
	std::string filename;

	double radius = 0.02;

	if (arg == 1)
	{
		std::cout << "Argument not founded " << std::endl;
		return 0;
	}

	filename = argv[1];

	if (arg == 3)
	{
		radius = atof(argv[2]) / 1000000;
	}

	std::ifstream input;
	input.open(filename.c_str());

	if (!input)
	{
		std::cout << "Cannot open data file" << std::endl;
		return 0;
	}

	PointType   point;
	std::size_t npoints;

	//Read GiD .post.msh format files NYI
	//     std::string GiD_Header;
	//     std::string Token;
	//     
	//     std::getline(input,GiD_Header)  //GiD Header
	//     std::getline(input,Token)       //CoCoordinates     

	input >> npoints;

	points = new PointType*[npoints];

	std::size_t pid;

	for (std::size_t i = 0; i < npoints; i++)
	{
		input >> pid;
		input >> point;

		points[i] = new PointType(point);
		points[i]->id = pid;
	}

	PointType min_point(*points[0]);
	PointType max_point(*points[0]);
	PointType mid_point;

	min_point.id = 0;
	max_point.id = 0;
	mid_point.id = 0;

	for (std::size_t i = 0; i < npoints; i++)
	{
		for (std::size_t j = 0; j < 3; j++)
		{
			if (min_point[j] >(*points[i])[j]) min_point[j] = (*points[i])[j];
			if (max_point[j] < (*points[i])[j]) max_point[j] = (*points[i])[j];
		}
	}

	for (std::size_t i = 0; i < Dim; i++)
	{
		mid_point.coord[i] = (max_point[i] + min_point[i]) / 2.00;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//  Output data Info                                                                                //
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	// Note: This version of main.cpp has a modified output format so as is easier to parse it.

	PointType*  search_point = &mid_point;

	std::size_t numsearch = 100000;
	std::size_t numsearch_nearest = numsearch * 100;

	std::cout << " min point : " << min_point << std::endl;
	std::cout << " max point : " << max_point << std::endl;
	std::cout << " search_point : " << *search_point << std::endl;
	std::cout << " search radius : " << radius << std::endl;
	std::cout << std::endl;

	std::cout << " Number of Points : " << npoints << std::endl;
	std::cout << " Number of Repetitions : " << numsearch << std::endl;
	std::cout << std::endl;

	std::cout << "SS\t\tGEN\tSIROMP\tSIRSER\tSNPOMP\tSNPSER\tNOFR\tNP" << std::endl;

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//  Data Setup                                                                                      //
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	PointType* allPoints;
	allPoints = new PointType[numsearch];

	std::size_t max_results = 100;// npoints;
	for (std::size_t i = 0; i < 1; i++)
	{
		allPoints[i] = *search_point;

		//radiusArray[i] = radius;

		//resultsArray[i] = new PtrPointType[npoints];
		//distancesArray[i] = new double[npoints];
	}


								   //Prepare the search point, search radius and resut arrays
	DistanceIterator distances = new double[npoints];
	PointIterator p_results = new PtrPointType[max_results];


	// Point-Based Search Structures
	std::vector<Point<3>> points_vector;
	for (std::size_t i = 0; i < npoints; i++)
	{
		points_vector.push_back(*(points[i]));
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//  Run tests                                                                                       //
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	double t0 = GetCurrentTime();
	PointsBins<Point<3>> bins(points_vector.begin(), points_vector.end());
	double t1 = GetCurrentTime();
	std::cout << "Points Bin" << "\t" << t1 - t0 << "\t";
	
	std::vector<PointsBins<Point<3>>::ResultType> results;
	PointsBins<Point<3>>::ResultType nearest_point_result;

	t0 = GetCurrentTime();
#pragma omp parallel for firstprivate(results)
	for (std::size_t i = 0; i < numsearch; i++) {
		results.clear();
		bins.SearchInRadius(*search_point, radius, results);
	}
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < numsearch; i++) {
		results.clear();
		bins.SearchInRadius(*search_point, radius, results);
	}
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << " \t";

	t0 = GetCurrentTime();

	for (std::size_t i = 0; i < numsearch_nearest; i++) {
		nearest_point_result = bins.SearchNearest(*search_point);
	}
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "  \t";

	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < numsearch_nearest; i++) {
		nearest_point_result = bins.SearchNearest(*search_point);
	}
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << " \t" << results.size() << "\t" << *nearest_point_result.Get();
	std::cout << std::endl;

	RunTestOMP<static_bins_type>("StaticBins", points, points + npoints, p_results, distances, max_results, allPoints, radius, numsearch, 1);
	RunTestOMP<dynamic_bins_type>("DynamicBins", points, points + npoints, p_results, distances, max_results, allPoints, radius, numsearch, 1);
	//RunTestOMP<kdtree_type>("KdTree\t", points, points + npoints, p_results, distances, max_results, allPoints, radius, numsearch, 10);
	//RunTestOMP<kdtree_average_split_type>("KdTreeAverage", points, points + npoints, resultsArray, distancesArray, max_results, allPoints, radiusArray, numsearch, 10);
	//RunTestOMP<kdtree_midpoint_split_type>("KdTreeMidpoint", points, points + npoints, resultsArray, distancesArray, max_results, allPoints, radiusArray, numsearch, 10);
	RunTestOMP<octree_type>("OcTree\t", points, points + npoints, p_results, distances, max_results, allPoints, radius, numsearch, 10);

	return 0;
}
