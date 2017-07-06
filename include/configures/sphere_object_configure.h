#pragma once

#include <limits>
#include <cmath>
#include <vector>

#include "point.h"
#include "sphere_object.h"

#include "contact_pair.h"

class SphereObjectConfigure {
public:

  static constexpr auto Epsilon   = std::numeric_limits<double>::epsilon();
  static constexpr auto Dimension = 3;

  typedef Point<3>          PointType;
  typedef SphereObject<3>   ObjectType;
  typedef SphereObject<3> * PointerType;

  typedef std::vector<ObjectType>             ObjectContainerType;
  typedef std::vector<PointerType>            ContainerType;

  typedef std::vector<PointerType>            ResultContainerType;
  typedef ContactPair<PointerType>            ContactPairType;

  typedef ContainerType::iterator             IteratorType;
  typedef ResultContainerType::iterator       ResultIteratorType;
  typedef std::vector<double>::iterator       DistanceIteratorType;

  typedef void * ContainerContactType;
  typedef void * IteratorContactType;

  /// Default consturctor
  SphereObjectConfigure(){};

  /// Default destructor
  virtual ~SphereObjectConfigure(){}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  /** Calculates the bounding box for the given object.
   * For this configuation file, the bounding box is the equal to the point given in 'rObject'.
   * @param rObject    Point for which the bounding box will be calculated.
   * @param rLowPoint  Lower point of the boundingbox.
   * @param rHighPoint Higher point of the boundingbox.
   */
  static inline void CalculateBoundingBox(const PointerType& rObject, Point<3>& rLowPoint, Point<3>& rHighPoint) {
    for(std::size_t d = 0; d < 3; d++) {
      rHighPoint[d] = rLowPoint[d] = (*rObject)[d];
    }
  }

  /** Calculates the bounding box for the given object extended with a Radius.
   * For this configuation file, the bounding box is the equal to the point given in 'rObject' + - a radius.
   * @param rObject    Point for which the bounding box will be calculated.
   * @param rLowPoint  Lower point of the boundingbox.
   * @param rHighPoint Higher point of the boundingbox.
   * @param Radius     The extension radius to be applied to the boundingbox.
   */
  static inline void CalculateBoundingBox(const PointerType& rObject, Point<3>& rLowPoint, Point<3>& rHighPoint, const double& Radius) {
    for(std::size_t d = 0; d < 3; d++) {
      rLowPoint[d] = (*rObject)[d] - Radius;
      rHighPoint[d] = (*rObject)[d] + Radius;
    }
  }

  /** Calculates the Center of the object.
   * @param rObject        Point for which the bounding box will be calculated.
   * @param rCentralPoint  The center point of the object.
   */
  static inline void CalculateCenter(const PointerType& rObject, Point<3>& rCentralPoint) {
    for(std::size_t d = 0; d < 3; d++) {
      rCentralPoint[d] = (*rObject)[d];
    }
  }

  /** Tests the intersection of two objects
   * For this configuation file, tests if the two points are the same within a Epsilon tolerance range.
   * @param  rObj_1 First point of the tests
   * @param  rObj_2 Second point of the tests
   * @return        Boolean indicating the result of the intersection test described.
   */
  static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2) {
    for(std::size_t i = 0; i < Dimension; i++) {
      if(std::fabs((*rObj_1)[i] - (*rObj_2)[i]) > Epsilon) {
        return false;
      }
    }

    return true;
  }

  /** Tests the intersection of two objects extended with a given radius.
   * For this configuation file, tests if the two points extended with a radius
   * are the same within a Epsilon tolerance range.
   * @param  rObj_1 First point of the tests
   * @param  rObj_2 Second point of the tests
   * @param  Radius The extension radius to be applied in the intersection.
   * @return        Boolean indicating the result of the intersection test described.
   */
  static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2, double Radius) {
    for(std::size_t i = 0; i < Dimension; i++) {
      if(std::fabs((*rObj_1)[i] - (*rObj_2)[i]) > Epsilon + Radius) {
        return false;
      }
    }

    return true;
  }

  /** Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
   * For this configuation file, tests if one point is inside the boundingbox
   * described by 'rLowPoint' and 'rHighPoint' within a Epsilon tolerance range.
   * @param  rObject    Point of the tests.
   * @param  rLowPoint  Lower point of the boundingbox.
   * @param  rHighPoint Higher point of the boundingbox.
   * @return            Boolean indicating the result of the intersection test described.
   */
  static inline bool IntersectionBox(const PointerType& rObject, const Point<3>& rLowPoint, const Point<3>& rHighPoint) {
    for(std::size_t i = 0; i < Dimension; i++) {
      if( (*rObject)[i] < rLowPoint[i] - Epsilon || (*rObject)[i] > rHighPoint[i] + Epsilon) {
        return false;
      }
    }

    return true;
  }

  /** Tests the intersection of one object with a boundingbox descrived by 'rLowPoint' and 'rHighPoint'.
   * For this configuation file, tests if one point extended by radius is inside the boundingbox
   * described by 'rLowPoint' and 'rHighPoint' within a Epsilon tolerance range.
   * @param  rObject    Point of the tests.
   * @param  rLowPoint  Lower point of the boundingbox.
   * @param  rHighPoint Higher point of the boundingbox.
   * @param  Radius     The extension radius to be applied in the intersection.
   * @return            Boolean indicating the result of the intersection test described.
   */
  static inline bool IntersectionBox(const PointerType& rObject, const Point<3>& rLowPoint, const Point<3>& rHighPoint, const double& Radius) {
    for(std::size_t i = 0; i < Dimension; i++) {
      if( ((*rObject)[i] + Radius) < rLowPoint[i] - Epsilon || ((*rObject)[i] - Radius) > rHighPoint[i] + Epsilon) {
        return false;
      }
    }

    return true;
  }

  /** Calculates the distance betwen two objects.
   * For this configuation file, calculates the euclidean distance between 'rObj_1' and 'rObj_2'.
   * # Performance
   * In C++11 'std::pow(T, int)' provides the optimal solution in terms of speed.
   * # References
   * (http://en.cppreference.com/w/cpp/numeric/math/pow)
   * (http://stackoverflow.com/questions/2940367)
   * @param rObj_1      First point.
   * @param rLowPoint   Lower point.
   * @param rHighPoint  Higher point of the boundingbox.
   * @param distance    The euclidean distance between 'rObj_1' and 'rObj_2'.
   */
  static inline void Distance(const PointerType& rObj_1, const PointerType& rObj_2, double& distance) {
    double pwdDistance = 0.0f;

    for(std::size_t i = 0; i < Dimension; i++) {
      pwdDistance += std::pow((*rObj_1)[i] - (*rObj_2)[i], 2);
    }

    distance = std::sqrt(pwdDistance);
  }

  /// Turns back information as a string.
  virtual std::string Info() const {
    return "Spatial Containers Configure for 'Points'";
  }

  /// Turns back data as a string.
  virtual std::string Data() const {
    return "Dimension: " + std::to_string(Dimension);
  }

  /// Prints object's information.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << Info() << std::endl;
  }

  /// Prints object's data.
  virtual void PrintData(std::ostream& rOStream) const {
    rOStream << Data() << Dimension << std::endl;
  }

protected:

private:

  /// Assignment operator.
  SphereObjectConfigure& operator=(SphereObjectConfigure const& rOther);

  /// Copy constructor.
  SphereObjectConfigure(SphereObjectConfigure const& rOther);
}; // Class SphereObjectConfigure

inline std::istream& operator >> (std::istream& rIStream, SphereObjectConfigure& rThis){
  return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const SphereObjectConfigure& rThis){
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
