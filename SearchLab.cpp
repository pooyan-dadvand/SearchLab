//System includes
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>

// Ugly fixes
#include <assert.h>
#define KRATOS_ERROR std::cout

// Kratos Independent
#define KRATOS_INDEPENDENT

// Interfaces
#include "interfaces/points_old_interface.h"
#include "interfaces/points_new_interface.h"
#include "interfaces/objects_old_interface.h"

constexpr std::size_t Dim = 3;

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

typedef Point<3> *        PtrPointType;
typedef PtrPointType *    PointVector;
typedef PtrPointType *    PointIterator;

typedef SphereObject<3> * PtrObjectType;
typedef PtrObjectType *   ObjectVector;
typedef PtrObjectType *   ObjectIterator;

typedef double* DistanceVector;
typedef double* DistanceIterator;

typedef Kratos::Bucket<Dim, Point<3>, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<Point<3>, Dim>>       BucketType;  //Bucket;
typedef Kratos::Bins<Dim, Point<3>, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<Point<3>, Dim>>         BinsStaticType;           //StaticBins;
typedef Kratos::BinsDynamic<Dim, Point<3>, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<Point<3>, Dim>>  BinsDynamicType; //DynamicBins;

int main(int arg, char* argv[]) {

	// Input data
	std::cout << std::setprecision(4) << std::fixed;

  //
	Point<3> ** points;
  SphereObject<3> ** objects;

	std::string filename;

	double radius = 0.02;

	if (arg == 1) {
		std::cout << "Argument not found" << std::endl;
		return 1;
	}

	filename = argv[1];

	if (arg == 3) {
		radius = atof(argv[2]) / 1000000;
	}

	std::ifstream input;
	input.open(filename.c_str());

	if (!input) {
		std::cout << "Cannot open data file" << std::endl;
		return 1;
	}

	Point<3> point;
  SphereObject<3> object;

	std::size_t npoints;

	input >> npoints;

	points = new Point<3>*[npoints];
  objects = new SphereObject<3>*[npoints];

	std::size_t pid;

	for (std::size_t i = 0; i < npoints; i++) {
		input >> pid;
		input >> point;

		points[i] = new Point<3>(point);
		points[i]->id = pid;
	}

	Point<3> min_point(*points[0]);
	Point<3> max_point(*points[0]);
	Point<3> mid_point;

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
	Point<3> & search_point = mid_point;

	std::size_t numsearch = 100000;
	std::size_t numsearch_nearest = numsearch * 100;

	std::cout << " min point : " << min_point << std::endl;
	std::cout << " max point : " << max_point << std::endl;
	std::cout << " search_point : " << search_point << std::endl;
	std::cout << " search radius : " << radius << std::endl;
	std::cout << std::endl;

	std::cout << " Number of Points : " << npoints << std::endl;
	std::cout << " Number of Repetitions : " << numsearch << std::endl;
	std::cout << std::endl;

	std::cout << "SS\t\tGEN\tSIROMP\tSIRSER\tSNPOMP\tSNPSER\tNOFR\tNP" << std::endl;

	// Data Setup
	Point<3> * allPoints = new Point<3>[numsearch];
  SphereObject<3> * allSpheres = new SphereObject<3>[numsearch];

	std::size_t max_results = 100;// npoints;
	for (std::size_t i = 0; i < 1; i++) {
		allPoints[i] = search_point;
	}

	//Prepare the search point, search radius and resut arrays
	DistanceIterator distances = new double[npoints];
	PointIterator p_results = new PtrPointType[max_results];
  ObjectIterator o_results = new PtrObjectType[max_results];

	// Point-Based Search Structures
	std::vector<Point<3>> points_vector;
	for (std::size_t i = 0; i < npoints; i++) {
		points_vector.push_back(*(points[i]));
	}

  // Point Interfaces
  // - New Interface
  PointsNew::RunTests<PointsBins<Point<3>>>("PointBins", points_vector, search_point, radius, numsearch, numsearch_nearest);

  // - Old Interface
	PointsOld::RunTests<BinsStaticType>("StaticBins", points, points + npoints, p_results, distances, max_results, allPoints, radius, numsearch, 1);
	PointsOld::RunTests<BinsDynamicType>("DynamicBins", points, points + npoints, p_results, distances, max_results, allPoints, radius, numsearch, 1);
  PointsOld::RunTests<OctreeType>("OcTree\t", points, points + npoints, p_results, distances, max_results, allPoints, radius, numsearch, 10);

  // Object Interfaces
  // - New Interface
  // TO BE FILLED

  // - Old Interface
  ObjectsOld::RunTests<BinsObjectStaticType>("BinsObjectStatic", objects, objects + npoints, o_results, distances, max_results, allPoints, radius, numsearch, 1);
  // RunTestsOldInterface<BinsObjectDynamicType>

	return 0;
}
