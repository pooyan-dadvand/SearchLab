#include <iostream>
#include <vector>

#include "timer.h"
#include "sphere_object.h"

#include "containers.h"
#include "configures/sphere_object_configure.h"
#include "spatial_containers/spatial_containers.h"

namespace ObjectsOld {

template<class SearchStructureType, class ObjectType, class IteratorType, class DistanceIterator>
static void RunTests(char const Title[], IteratorType PBegin, IteratorType PEnd, IteratorType results0, DistanceIterator distances0, std::size_t MaxResults, ObjectType* allObjects, double radius0, std::size_t numsearch, std::size_t bucket_size) {
	double t0, t1;
	std::size_t n;

	std::size_t npoints = PEnd - PBegin;
	std::vector<std::size_t> numresArray(numsearch);
	Entities::PtrObjectType objectToSearch = &allObjects[0];

	std::size_t numsearch_nearest = 100 * numsearch;

	t0 = GetCurrentTime();
  Kratos::BinsObjectDynamic<SphereObjectConfigure> nodes_tree(PBegin, PEnd);
	t1 = GetCurrentTime();
	std::cout << Title << "\t" << t1 - t0 << "\t";

	t0 = GetCurrentTime();
	#pragma omp parallel for
	for (std::size_t i = 0; i < 1; i++) {
		n = nodes_tree.SearchObjectsInRadius(objectToSearch, radius0, results0, MaxResults);
  }
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	t0 = GetCurrentTime();
	for (std::size_t i = 0; i < 1; i++) {
		n = nodes_tree.SearchObjectsInRadius(objectToSearch, radius0, results0, distances0, MaxResults);
  }
	t1 = GetCurrentTime();
	std::cout << t1 - t0 << "\t";

	ObjectType** PNearestArray = new ObjectType*[numsearch];
	std::vector<double> distancesArrayNearest(numsearch);
	ObjectType* PNearest = PNearestArray[0];

	// t0 = GetCurrentTime();
	// for (std::size_t i = 0; i < 100; i++) {
	// 	nodes_tree.SearchNearestPoint(allPoints, numsearch, PNearestArray, distancesArrayNearest);
  // }
	// t1 = GetCurrentTime();

	std::cout << "NaN" << "\t";
	double distance = 0;


	// t0 = GetCurrentTime();
	// for (std::size_t i = 0; i < numsearch_nearest; i++) {
	// 	PNearest = nodes_tree.SearchNearestPoint(point0, distance);
  // }
	// t1 = GetCurrentTime();
	std::cout << "NaN" << "\t";

	std::cout << n << "\t" << "Does not apply" << std::endl;
};

} // Namespace ObjectsOld
