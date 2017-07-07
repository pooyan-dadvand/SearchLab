#pragma once

#include "point.h"
#include "sphere_object.h"

#include "configures/sphere_object_configure.h"
#include "spatial_containers/spatial_containers.h"

constexpr std::size_t Dim = 3;

// Auxiliar function. Not very sure why this is here.
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

// Entities used in the tests
namespace Entities {
  typedef Point<3>          PointType;
  typedef Point<3> *        PtrPointType;
  typedef PtrPointType *    PointIterator;

  typedef SphereObject<3>   ObjectType;
  typedef SphereObject<3> * PtrObjectType;
  typedef std::vector<ObjectType>::iterator   ObjectIterator;

  typedef std::vector<double>::iterator DistanceIterator;
}

// Containers used in the tests
namespace Containers {

  typedef Entities::PtrPointType * PointVector;
  typedef Entities::PtrObjectType * ObjectVector;
  typedef double* DistanceVector;

  // Bucket ( Can't see any use besides the octree )
  typedef Kratos::Bucket<
    Dim,
    Entities::PointType,
    Containers::PointVector,
    Entities::PtrPointType,
    Entities::PointIterator,
    Entities::DistanceIterator,
    PointDistance2<
      Entities::PointType,
      Dim >
  > BucketType;

  // Static Bins
  typedef Kratos::Bins<
    Dim,
    Entities::PointType,
    Containers::PointVector,
    Entities::PtrPointType,
    Entities::PointIterator,
    Entities::DistanceIterator,
    PointDistance2<
      Entities::PointType,
      Dim >
  > BinsStaticType;

  // Dynamic Bins
  typedef Kratos::BinsDynamic<
    Dim,
    Entities::PointType,
    Containers::PointVector,
    Entities::PtrPointType,
    Entities::PointIterator,
    Entities::DistanceIterator,
    PointDistance2<
      Entities::PointType,
      Dim >
  > BinsDynamicType;

  // Octree Bins (Static);
  typedef Kratos::Tree<
    Kratos::OCTreePartition<BinsStaticType>
  > BinsStaticOctreeType;

  // Octree Bins (Dynamic);
  typedef Kratos::Tree<
    Kratos::OCTreePartition<BinsDynamicType>
  > BinsDynamicOctreeType;

  // Static Objects Bins
  typedef Kratos::BinsObjectStatic<
    SphereObjectConfigure
  > BinsObjectStaticType;

  // Dynamic Objects Bins
  typedef Kratos::BinsObjectDynamic<
    SphereObjectConfigure
  > BinsObjectDynamicType;
}
