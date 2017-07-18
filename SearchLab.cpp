//System includes
#include <fstream>
#include <iostream>
#include <ctime>
#include <string>
#include <cstdlib>
#include <iomanip>

#if __cplusplus < 201103L
#error message( "This library needs a C++11 compiler")
#endif

// Point
#include "point.h"

// Ugly fixes
#include <assert.h>

// Kratos Independent
#define KRATOS_INDEPENDENT

#ifndef KRATOS_INDEPENDENT
#define KRATOS_ERROR std::cout

// Kratos includes
#include "spatial_containers/spatial_containers.h"
#endif // KRATOS_INDEPENDENT

#include "points_bins.h"
#include "points_hash.h"

double GetCurrentTime() {
#ifndef _OPENMP
	return std::clock() / static_cast<double>(CLOCKS_PER_SEC);
#else
	return omp_get_wtime();
#endif
}

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

template<class SearchStructureType, class PointType, class IteratorType, class DistanceIterator>
void RunTestsOldInterface(char const Title[], IteratorType PBegin, IteratorType PEnd, IteratorType results0, DistanceIterator distances0, std::size_t MaxResults, PointType* allPoints, double radius0, std::size_t numsearch, std::size_t bucket_size) {
	double t0, t1;
	std::size_t n;

	std::size_t npoints = PEnd - PBegin;
	std::vector<std::size_t> numresArray(numsearch);
	PointType& point0 = allPoints[0];

	std::size_t numsearch_nearest = 10 * numsearch;

	t0 = GetCurrentTime();
	SearchStructureType nodes_tree(PBegin, PEnd);
	t1 = GetCurrentTime();
	std::cout << Title << "\t" << t1 - t0 << "\t";

	t0 = GetCurrentTime();
  #pragma omp parallel for
	for (std::size_t i = 0; i < numsearch; i++) {
		n = nodes_tree.SearchInRadius(point0, radius0, results0, distances0, MaxResults);
  }
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < numsearch; i++) {
		n = nodes_tree.SearchInRadius(point0, radius0, results0, distances0, MaxResults);
  }
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	PointType** PNearestArray = new PointType*[numsearch];
	std::vector<double> distancesArrayNearest(numsearch);
	PointType* PNearest = PNearestArray[0];

	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < 10; i++) {
		nodes_tree.SearchNearestPoint(allPoints, numsearch, PNearestArray, distancesArrayNearest);
  }
	t1 = GetCurrentTime();

	std::cout << t1 - t0 << "\t";
	double distance = 0;


	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < numsearch_nearest; i++) {
		PNearest = nodes_tree.SearchNearestPoint(point0, distance);
  }
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	std::cout << n << "\t" << *PNearestArray[0] << std::endl;
};

template<class BinsType>
void RunTestsNewInterface(char const Title[], std::vector<Point<3>> & points_vector, const Point<3> & search_point, double radius, std::size_t numsearch, std::size_t numsearch_nearest) {
  double t0 = GetCurrentTime();
  BinsType bins(points_vector.begin(), points_vector.end());
  double t1 = GetCurrentTime();
  std::cout << "Points Bin" << "\t" << t1 - t0 << "\t";

  std::vector<PointsBins<Point<3>>::ResultType> results;
  PointsBins<Point<3>>::ResultType nearest_point_result;

  t0 = GetCurrentTime();
  #pragma omp parallel for firstprivate(results)
  for (std::size_t i = 0; i < numsearch; i++) {
    results.clear();
    bins.SearchInRadius(search_point, radius, results);
  }
  t1 = GetCurrentTime();
  std::cout << t1 - t0 << "\t";

  t0 = GetCurrentTime();
  for (std::size_t i = 0; i < numsearch; i++) {
    results.clear();
    bins.SearchInRadius(search_point, radius, results);
  }
  t1 = GetCurrentTime();
  std::cout << t1 - t0 << "\t";

  t0 = GetCurrentTime();

  #pragma omp parallel for firstprivate(nearest_point_result)
  for (std::size_t i = 0; i < numsearch_nearest; i++) {
    nearest_point_result = bins.SearchNearest(search_point);
  }
  t1 = GetCurrentTime();
  std::cout << t1 - t0 << "\t";

  t0 = GetCurrentTime();
  for (std::size_t i = 0; i < numsearch_nearest; i++) {
    nearest_point_result = bins.SearchNearest(search_point);
  }
  t1 = GetCurrentTime();
  std::cout << t1 - t0 << "\t" << results.size() << "\t" << *nearest_point_result.Get();
  std::cout << std::endl;
}

int RunPointSearchComparison(std::string Filename, double Radius) {
	constexpr std::size_t Dim = 3;

	typedef Point<Dim>       PointType;

	typedef PointType *      PtrPointType;
	typedef PtrPointType *   PointVector;
	typedef PtrPointType *   PointIterator;

	typedef double* DistanceVector;
	typedef double* DistanceIterator;

#ifndef KRATOS_INDEPENDENT
	typedef Kratos::Bucket<Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<PointType, Dim>>       bucket_type;  //Bucket;
	typedef Kratos::Bins<Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<PointType, Dim>>         StaticBinsType;           //StaticBins;
	typedef Kratos::BinsDynamic<Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<PointType, Dim>>  DynamicBinsType; //DynamicBins;

	typedef Kratos::Tree< Kratos::KDTreePartition<bucket_type>>                  kdtree_type;                //Kdtree;
	typedef Kratos::Tree< Kratos::KDTreePartitionAverageSplit<bucket_type>>      kdtree_average_split_type;  //Kdtree;
	typedef Kratos::Tree< Kratos::KDTreePartitionMidPointSplit<bucket_type>>     kdtree_midpoint_split_type; //Kdtree;
	typedef Kratos::Tree< Kratos::OCTreePartition<bucket_type>>                  OctreeType;                 //Octree;
	typedef Kratos::Tree< Kratos::KDTreePartitionMidPointSplit<StaticBinsType>>  kdtree_StaticBinsType;      //KdtreeBins;
	typedef Kratos::Tree< Kratos::OCTreePartition<StaticBinsType>>               octree_StaticBinsType;      //OctreeBins;
	typedef Kratos::Tree< Kratos::KDTreePartitionMidPointSplit<DynamicBinsType>> kdtree_DynamicBinsType;     //KdtreeBins;
	typedef Kratos::Tree< Kratos::OCTreePartition<DynamicBinsType>>              octree_bins_type;           //OctreeBins;
#endif // KRATOS_INDEPENDENT

	// Input data
	std::cout << std::setprecision(4) << std::fixed;

	PointVector points;

	std::ifstream input;
	input.open(Filename.c_str());

	if (!input) {
		std::cout << "Cannot open data file" << std::endl;
		return 0;
	}
	std::cout << "Comparison for " << Filename << std::endl;
	PointType   point;
	std::size_t npoints;

	input >> npoints;
	points = new PointType*[npoints];

	std::size_t pid;

	for (std::size_t i = 0; i < npoints; i++) {
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

	for (std::size_t i = 0; i < npoints; i++) {
		for (std::size_t j = 0; j < 3; j++) {
			if (min_point[j] > (*points[i])[j]) min_point[j] = (*points[i])[j];
			if (max_point[j] < (*points[i])[j]) max_point[j] = (*points[i])[j];
		}
	}

	for (std::size_t i = 0; i < Dim; i++) {
		mid_point.coord[i] = (max_point[i] + min_point[i]) / 2.00;
	}

	// Output data Info
	PointType & search_point = mid_point;

	std::size_t numsearch = 100000;
	std::size_t numsearch_nearest = numsearch * 10;

	std::cout << " min point : " << min_point << std::endl;
	std::cout << " max point : " << max_point << std::endl;
	std::cout << " search_point : " << search_point << std::endl;
	std::cout << " search radius : " << Radius << std::endl;
	std::cout << std::endl;

	std::cout << " Number of Points : " << npoints << std::endl;
	std::cout << " Number of Repetitions : " << numsearch << std::endl;
	std::cout << std::endl;

	std::cout << "SS\t\tGEN\tSIROMP\tSIRSER\tSNPOMP\tSNPSER\tNOFR\tNP" << std::endl;

	// Data Setup
	PointType * allPoints;
	allPoints = new PointType[numsearch];

	std::size_t max_results = npoints;
	for (std::size_t i = 0; i < 1; i++) {
		allPoints[i] = search_point;
	}

	//Prepare the search point, search radius and resut arrays
	DistanceIterator distances = new double[npoints];
	PointIterator p_results = new PtrPointType[max_results];

	// Point-Based Search Structures
	std::vector<Point<3>> points_vector;
	for (std::size_t i = 0; i < npoints; i++) {
		points_vector.push_back(*(points[i]));
	}

  // New Interface
	RunTestsNewInterface<PointsBins<Point<3>>>("PointBins", points_vector, search_point, Radius, numsearch, numsearch_nearest);

  // Old Interface
#ifndef KRATOS_INDEPENDENT
	RunTestsOldInterface<StaticBinsType>("StaticBins", points, points + npoints, p_results, distances, max_results, allPoints, Radius, numsearch, 1);
	RunTestsOldInterface<DynamicBinsType>("DynamicBins", points, points + npoints, p_results, distances, max_results, allPoints, Radius, numsearch, 1);
  RunTestsOldInterface<OctreeType>("OcTree\t", points, points + npoints, p_results, distances, max_results, allPoints, Radius, numsearch, 10);

  //RunTestsOldInterface<kdtree_type>("KdTree\t", points, points + npoints, p_results, distances, max_results, allPoints, radius, numsearch, 10);
	//RunTestsOldInterface<kdtree_average_split_type>("KdTreeAverage", points, points + npoints, resultsArray, distancesArray, max_results, allPoints, radiusArray, numsearch, 10);
	//RunTestsOldInterface<kdtree_midpoint_split_type>("KdTreeMidpoint", points, points + npoints, resultsArray, distancesArray, max_results, allPoints, radiusArray, numsearch, 10);
#endif // KRATOS_INDEPENDENT

	return 0;
}

int main(int arg, char* argv[]) {
	std::string filename;

	double radius = 0.01;


	if (arg > 1) {
		std::cout << "Argument not founded " << std::endl;
		filename = argv[1];

		if (arg == 3) {
			radius = atof(argv[2]) / 1000000;
		}


		return 0;
	}

	filename = "genericCube100x100x100.5051.pts";
	RunPointSearchComparison(filename, radius);
	filename = "offsetCube79x79x79.1603.pts";
	RunPointSearchComparison(filename, radius);
	filename = "clusterCube6x6x6X4913.490.pts";
	RunPointSearchComparison(filename, radius);
	filename = "line100000.5.pts";
	RunPointSearchComparison(filename, radius);

	return 0;
}
