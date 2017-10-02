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

int RunPointSearchComparison(std::string Filename, double Radius, int sizeChunksX, int sizeChunksY, int sizeChunksZ) {

	int mpi_rank, mpi_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	// Input data
	std::cout << std::setprecision(4) << std::fixed;

	Point<3> ** points;
  Point<3> point;
  SphereObject<3> object;

	std::ifstream input;
	input.open(Filename.c_str());

	if (!input) {
		if(!mpi_rank) {
			std::cout << "Cannot open data file" << std::endl;
		}
		return 1;
	}

	if(!mpi_rank) {
		std::cout << "Comparison for " << Filename << std::endl;
	}

	std::size_t npoints;
	input >> npoints;

	int sizeBlock = 2;									// Num of blocks per chunk
	std::size_t npointstot = npoints;

	points = new Point<3>*[npointstot];
  // std::vector<Entities::PtrObjectType> objects(npointstot);

	std::size_t pid;

	int Sdx = 4;
	int Sdy = 4;
	int Sdz = 3;

	int Sdi = 0;

	for(int i = 48; i < mpi_size; i*=2) {
		if(Sdi == 0) Sdx *= 2;
		if(Sdi == 1) Sdy *= 2;
		if(Sdi == 2) Sdz *= 2;
		Sdi = (Sdi + 1) % 3;
	}

	int Idx = mpi_rank % Sdx;
	int Idy = (mpi_rank % (Sdx * Sdy)) / Sdx;
	int Idz = mpi_rank / (Sdx * Sdy);

	// std::cout << mpi_rank << " " << Idx << " " << Idy << " " << Idz << std::endl;

	for(std::size_t i = 0; i < npoints; i++) {

		input >> pid;
		input >> point;

    for(std::size_t d = 0; d < 3; d++) {
      object[d] = point[d];
    }

    object.radius = 0.5/npoints;

		// for(int z = 0; z < sizeChunksZ; z++) {
		// 	for(int y = 0; y < sizeChunksY; y++) {
		// 		for(int x = 0; x < sizeChunksX; x++) {
		// 			int index = i + x * npoints + y * sx * npoints + z * sx * sy * npoints;

		int index = i;
		points[index] = new Point<3>(point);
		points[index]->id = pid;
		(*points[index])[0] += Idx;
		(*points[index])[1] += Idy;
		(*points[index])[2] += Idz;

		// objects[index] = new SphereObject<3>(object);
		// objects[index]->id = pid;
		// objects[index]->radius = 0.5/npoints;
		// 		}
		// 	}
		// }
	}

	Point<3> min_point(*points[0]);
	Point<3> max_point(*points[0]);
	Point<3> mid_point;
  SphereObject<3> mid_object;

	// min_point.id = 0;
	// max_point.id = 0;
	// mid_point.id = 0;

	for (std::size_t i = 0; i < npointstot; i++) {
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

	std::size_t numsearch = npointstot;
	std::size_t numsearch_nearest = npointstot * 10;

	if(!mpi_rank) {
		std::cout << " min point : " << min_point << std::endl;
		std::cout << " max point : " << max_point << std::endl;
		std::cout << " search_point : " << search_point << std::endl;
		std::cout << " search radius : " << Radius << std::endl;
		std::cout << std::endl;

		std::cout << " Number of Points : " << npoints << std::endl;
		std::cout << " Number of Repetitions : " << numsearch << std::endl;
		std::cout << std::endl;

		std::cout << "SS\t\tGEN\tSIROMP\tSIRSER\tSNPOMP\tSNPSER\tNOFR\tNP" << std::endl;
	}

	// // Data Setup
	// Point<3> * allPoints = new Point<3>[numsearch];
  // SphereObject<3> * allSpheres = new SphereObject<3>[numsearch];

	// std::size_t max_results = npoints;
	// for (std::size_t i = 0; i < 1; i++) {
	// 	allPoints[i] = search_point;
  //   allSpheres[i] = search_object;
	// }

	//Prepare the search point, search radius and resut arrays

#ifdef USE_KRATOS
  std::vector<Entities::PtrObjectType> objectResults(max_results);
  std::vector<double> resultDistances(max_results);


	double * distances = new double[npoints];
	Entities::PointIterator p_results = new Entities::PtrPointType[max_results];
#endif

	// Point-Based Search Structures
	std::vector<Point<3>> points_vector;
	for (std::size_t i = 0; i < npointstot; i++) {
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

	MPI_Init(&arg, &argv);

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

	// filename = "../cases/genericCube3x3x3.250000.pts";
	filename = "../cases/randomCube2000000.pts";
	RunPointSearchComparison(filename, radius, 4, 4, 3);
	
	// filename = "../cases/offsetCube79x79x79.1603.pts";
	// RunPointSearchComparison(filename, radius);
	// filename = "../cases/clusterCube6x6x6X4913.490.pts";
	// RunPointSearchComparison(filename, radius);
	// filename = "../cases/line100000.5.pts";
	// RunPointSearchComparison(filename, radius);

	MPI_Finalize();

	return 0;
}
