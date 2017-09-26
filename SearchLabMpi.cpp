//System includes
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>

// Kratos
#ifdef USE_KRATOS
	// Ugly fixes
	#include <assert.h>
	#define KRATOS_ERROR std::cout

	#include "interfaces/points_old_interface.h"
	#include "interfaces/objects_old_interface.h"
#endif

// Containers
#include "containers.h"
#include "parallel_bins.h"

// Interfaces
#include "interfaces/points_new_interface_mpi.h"

//Octree includes
#include "custom_utilities/octree_driver.h"
#include "custom_utilities/octree_binary_cell.h"

int RunPointSearchComparison(std::string Filename, double Radius) {

	// Input data
	std::cout << std::setprecision(4) << std::fixed;

	Point<3> ** points;
  Point<3> point;
  SphereObject<3> object;

	std::ifstream input;
	input.open(Filename.c_str());

	if (!input) {
		std::cout << "Cannot open data file" << std::endl;
		return 1;
	}

	std::cout << "Comparison for " << Filename << std::endl;

	std::size_t npoints;

	input >> npoints;

	points = new Point<3>*[npoints];
  std::vector<Entities::PtrObjectType> objects(npoints);

	std::size_t pid;

	for(std::size_t i = 0; i < npoints; i++) {
		input >> pid;
		input >> point;

    for(std::size_t d = 0; d < 3; d++) {
      object[d] = point[d];
    }
    object.radius = 0.5/npoints;

		points[i] = new Point<3>(point);
		points[i]->id = pid;

    objects[i] = new SphereObject<3>(object);
    objects[i]->id = pid;
    objects[i]->radius = 0.5/npoints;
	}

	Point<3> min_point(*points[0]);
	Point<3> max_point(*points[0]);
	Point<3> mid_point;
  SphereObject<3> mid_object;

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
    mid_object.coord[i] = (max_point[i] + min_point[i]) / 2.00;
	}

  mid_object.radius = 0.5/npoints;

	// Output data Info
	Point<3> & search_point = mid_point;
  SphereObject<3> & search_object = mid_object;

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
	Point<3> * allPoints = new Point<3>[numsearch];
  SphereObject<3> * allSpheres = new SphereObject<3>[numsearch];

	std::size_t max_results = npoints;
	for (std::size_t i = 0; i < 1; i++) {
		allPoints[i] = search_point;
    allSpheres[i] = search_object;
	}

	//Prepare the search point, search radius and resut arrays

  std::vector<Entities::PtrObjectType> objectResults(max_results);
  std::vector<double> resultDistances(max_results);

	double * distances = new double[npoints];
	Entities::PointIterator p_results = new Entities::PtrPointType[max_results];

	// Point-Based Search Structures
	std::vector<Point<3>> points_vector;
	for (std::size_t i = 0; i < npoints; i++) {
		points_vector.push_back(*(points[i]));
	}

  // Point Interfaces
  // - New Interface
  PointsNew::RunTests<ParallelBins<PointsBins<Point<3>>>>("PointBins", points_vector, search_point, Radius, numsearch, numsearch_nearest);

  // - Old Interface
#ifdef USE_KRATOS
	PointsOld::RunTests<Containers::BinsStaticType>("StaticBins", points, points + npoints, p_results, resultDistances.begin(), max_results, allPoints, Radius, numsearch, 1);
	PointsOld::RunTests<Containers::BinsDynamicType>("DynamicBins", points, points + npoints, p_results, resultDistances.begin(), max_results, allPoints, Radius, numsearch, 1);
  PointsOld::RunTests<Containers::OctreeType>("OctTree\t", points, points + npoints, p_results, resultDistances.begin(), max_results, allPoints, Radius, numsearch, 10);
  // PointsOld::RunTests<Containers::BinsDynamicOctreeType>("OcTreeDynamic\t", points, points + npoints, p_results, resultDistances.begin(), max_results, allPoints, Radius, numsearch, 10);
#endif

  // Object Interfaces
  // - New Interface
  // TO BE FILLED

  // - Old Interface
#ifdef USE_KRATOS
  ObjectsOld::RunTests<Containers::BinsObjectStaticType>("StaticObjects", objects.begin(), objects.end(), objectResults.begin(), resultDistances.begin(), max_results, allSpheres, Radius, numsearch, 1);
  ObjectsOld::RunTests<Containers::BinsObjectDynamicType>("DynamicObjects", objects.begin(), objects.end(), objectResults.begin(), resultDistances.begin(), max_results, allSpheres, Radius, numsearch, 1);
#endif

	return 0;
}

int main(int arg, char* argv[]) {

  double center[ 3 ] = {0.5,0.5,0.5};
	double radius = 1.5;
  int level = 12;
  

  //TESTS WITH OCTREE DRIVER
	std::string filename;
	  if (arg > 1) {
		  std::cout << "Argument not founded " << std::endl;
		  filename = argv[1];
		  if (arg == 3) {
  	    radius = atof(argv[2]) / 1000000;
		  }
    return 0;
  }
  //filename = "../cases/genericCube2x2x2.500000.pts";
	//filename = "../cases/genericCube10x10x10.55556.pts";
	//filename = "../cases/genericCube100x100x100.5051.pts";
	//filename = "../cases/line100000.5.pts";
	filename = "../cases/offsetCube79x79x79.1603.pts";
  RunPointSearchOctree( filename , radius , center , level , arg , argv );

	return 0;
}











