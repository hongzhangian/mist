// Copyright (c) 2014, Yan, Hongzhang. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#ifndef MIST_MATH_MATH_UTIL_H_
#define MIST_MATH_MATH_UTIL_H_

#include <vector>
#include <list>
#include <utility>
#include <algorithm>
#include <limits>

namespace mist {
namespace math {

class Radian;
class Degree;
class Math;

class Angle {
 public:
  explicit Angle(Real angle) : angle_(angle) {}
  operator Radian() const {
    return Radian(Math::AngleUnitsToRadians(angle_));
  }
  operator Degree() const {
    return Degree(Math::AngleUnitsToDegrees(angle_));
  }

 private:
  Real angle_;
};

class Math {
 public:
  // The angular units used by the API. This functionality is now deprecated in
  // favor of discreet angular unit types ( see Degree and Radian above ).
  // The only place this functionality is actually still used is when parsing
  // files. Search for usage of the Angle class for those instances
  enum AngleUnit {
    AU_DEGREE,
    AU_RADIAN
  };

  // Default constructor.
  // trig_table_size Optional parameter to set the size of the tables used to
  // implement Sin, Cos, Tan
  explicit Math(unsigned int trig_table_size = 4096);
  ~Math();

  static inline int IAbs(int iValue) {
    return ( iValue >= 0 ? iValue : -iValue );
  }
  static inline int ICeil(float value) {
    return  static_cast<int>(ceil(value));
  }
  static inline int IFloor(float value) {
    return static_cast<int>(floor(value));
  }
  static int ISign(int iValue);

  static inline Real Abs(Real value) {
    return Real(fabs(value));
  }
  static inline Degree Abs(const Degree& dValue) {
    return Degree(fabs(dValue.valueDegrees()));
  }
  static inline Radian Abs(const Radian& rValue) {
    return Radian(fabs(rValue.ValueAsRadians()));
  }

  static inline Real Ceil(Real value) { return Real(ceil(value)); }
  static inline Real Floor(Real value) { return Real(floor(value)); }

  static inline bool IsNaN(Real f) {
    // std::isnan() is C99, not supported by all compilers
    // However NaN always fails this next test, no other number does.
    return f != f;
  }

  static inline Real Exp(Real value) { return Real(exp(value)); }

  static inline Real Log(Real value) { return Real(log(value)); }
  // Stored value of log(2) for frequent use
  static const Real LOG2;
  static inline Real Log2(Real value) { return Real(log(value)/LOG2); }
  static inline Real LogN(Real base, Real value) {
    return Real(log(value)/log(base));
  }
  static inline Real Pow(Real base, Real exponent) {
    return Real(pow(base, exponent));
  }

  static Real Sign(Real value);
  static inline Radian Sign(const Radian& radian) {
    return Radian(Sign(radian.ValueAsRadians()));
  }
  static inline Degree Sign(const Degree& degree) {
    return Degree(Sign(degree.ValueAsDegrees()));
  }

  static inline Real Sin(const Radian& value, bool use_tables = false) {
    return (!use_tables) ?
           Real(sin(value.ValueAsRadians())) :
           SinTable(value.ValueAsRadians());
  }
  static inline Real Sin(Real value, bool use_tables = false) {
    return (!use_tables) ? Real(sin(value)) : SinTable(value);
  }
  static inline Real Cos(const Radian& value, bool use_tables = false) {
    return (!use_tables) ?
           Real(cos(value.ValueAsRadians())) :
           SinTable(value.ValueAsRadians() + HALF_PI);
  }

  static inline Real Cos(Real value, bool use_tables = false) {
    return (!use_tables) ? Real(cos(value)) : SinTable(value + HALF_PI);
  }
  static inline Real Tan(const Radian& value, bool use_tables = false) {
    return (!use_tables) ?
           Real(tan(value.ValueAsRadians())) :
           TanTable(value.ValueAsRadians());
  }
  static inline Real Tan(Real value, bool use_tables = false) {
    return (!use_tables) ? Real(tan(value)) : TanTable(value);
  }

  static Radian ACos(Real value);
  static Radian ASin(Real value);
  static inline Radian ATan(Real value) { return Radian(atan(value)); }
  static inline Radian ATan2(Real x, Real y) { return Radian(atan2(x, y)); }

  static inline Real Sqr(Real value) { return value*value; }
  static inline Real Sqrt(Real value) { return Real(sqrt(value)); }
  static inline Radian Sqrt(const Radian& value) {
    return Radian(sqrt(value.ValueAsRadians()));
  }
  static inline Degree Sqrt(const Degree& value) {
    return Degree(sqrt(value.ValueAsDegrees()));
  }
  // Inverse square root i.e. 1 / Sqrt(x), good for vector
  static Real InvSqrt(Real value);

  // Generate a random number of unit length, in the range from [0,1].
  static Real UnitRandom();
  // Generate a random number within the range provided.
  static Real RangeRandom(Real low, Real high);
  // Generate a random number in the range [-1,1].
  static Real SymmetricRandom();


  static inline Real DegreesToRadians(Real degrees) {
    return degrees * deg_to_rad_;
  }
  static inline Real RadiansToDegrees(Real radians) {
    return radians * rad_to_deg_;
  }

  // These functions used to set the assumed angle units (radians or degrees)
  // expected when using the Angle type. You can set this directly after
  // creating a new Root, and also before/after resource creation, depending on
  // whether you want the change to affect resource files.
  static void SetAngleUnit(AngleUnit unit);
  // Get the unit being used for angles.
  static AngleUnit GetAngleUnit(void);

  // Convert from the current AngleUnit to radians.
  static Real AngleUnitsToRadians(Real units);
  // Convert from radians to the current AngleUnit .
  static Real RadiansToAngleUnits(Real radians);
  // Convert from the current AngleUnit to degrees.
  static Real AngleUnitsToDegrees(Real units);
  // Convert from degrees to the current AngleUnit.
  static Real DegreesToAngleUnits(Real degrees);

  // Checks whether a given point is inside a triangle, in a
  // The vertices of the triangle must be given in either
  // trigonometrical (anticlockwise) or inverse trigonometrical
  // (clockwise) order.
  static bool PointInTri2D(const Vector2& p,
                           const Vector2& a,
                           const Vector2& b,
                           const Vector2& c);

  // Checks whether a given 3D point is inside a triangle.
  // The vertices of the triangle must be given in either
  // trigonometrical (anticlockwise) or inverse trigonometrical
  // (clockwise) order, and the point must be guaranteed to be in the
  // normal The triangle plane's normal (passed in rather than calculated
  // on demand since the caller may already have it)

  static bool PointInTri3D(const Vector3& p,
                           const Vector3& a,
                           const Vector3& b,
                           const Vector3& c,
                           const Vector3& normal);

  // Ray / plane intersection, returns boolean result and distance.
  static std::pair<bool, Real> Intersects(const Ray& ray, const Plane& plane);

  // Ray / sphere intersection, returns boolean result and distance.
  static std::pair<bool, Real> Intersects(const Ray& ray,
                                          const Sphere& sphere,
                                          bool discard_inside = true);

  // Ray / box intersection, returns boolean result and distance.
  static std::pair<bool, Real> Intersects(const Ray& ray,
                                          const AxisAlignedBox& box);

  // Ray / box intersection, returns boolean result and two intersection
  // distance. d1, d2 are real pointer to retrieve the near intersection
  // distance from the ray origin, maybe NULL which means don't care
  // Return if the ray is intersects the box, true is returned, and the near
  // intersection distance is return by d1, the far intersection distance is
  // return by d2. Guarantee 0 <= d1 <= d2.
  static bool Intersects(const Ray& ray,
                         const AxisAlignedBox& box,
                         Real* d1,
                         Real* d2);

  // Ray / triangle intersection, returns boolean result and distance.
  static std::pair<bool, Real> Intersects(const Ray& ray,
                                          const Vector3& a,
                                          const Vector3& b,
                                          const Vector3& c,
                                          const Vector3& normal,
                                          bool positive_side = true,
                                          bool negative_side = true);

  // Ray / triangle intersection, returns boolean result and distance.
  static std::pair<bool, Real> Intersects(const Ray& ray,
                                          const Vector3& a,
                                          const Vector3& b,
                                          const Vector3& c,
                                          bool positive_side = true,
                                          bool negative_side = true);

  // Sphere / box intersection test.
  static bool Intersects(const Sphere& sphere, const AxisAlignedBox& box);

  // Plane / box intersection test.
  static bool Intersects(const Plane& plane, const AxisAlignedBox& box);

  // * Ray / convex plane list intersection test.
  static std::pair<bool, Real> Intersects(const Ray& ray,
                                          const vector<Plane>::type& planeList,
                                          bool normal_is_outside);

  // Ray / convex plane list intersection test.
  static std::pair<bool, Real> Intersects(const Ray& ray,
                                          const list<Plane>::type& planeList,
                                          bool normal_is_outside);

  // Sphere / plane intersection test.
  static bool Intersects(const Sphere& sphere, const Plane& plane);

  // Compare 2 reals, using tolerance for inaccuracies.
  static bool RealEqual(Real a, Real b,
                        Real tolerance = std::numeric_limits<Real>::epsilon());

  // Calculates the tangent space vector for a given set of positions / texture
  // coords.
  static Vector3 CalculateTangentSpaceVector(const Vector3& position1,
                                             const Vector3& position2,
                                             const Vector3& position3,
                                             Real u1, Real v1,
                                             Real u2, Real v2,
                                             Real u3, Real v3);

  // Build a reflection matrix for the passed in plane.
  static Matrix4 BuildReflectionMatrix(const Plane& p);
  // Calculate a face normal, including the w component which is the offset from
  // the origin.
  static Vector4 CalculateFaceNormal(const Vector3& v1,
                                     const Vector3& v2,
                                     const Vector3& v3);
  // Calculate a face normal, no w-information.
  static Vector3 CalculateBasicFaceNormal(const Vector3& v1,
                                          const Vector3& v2,
                                          const Vector3& v3);
  // Calculate a face normal without normalize, including the w component
  // which is the offset from the origin.
  static Vector4 CalculateFaceNormalWithoutNormalize(const Vector3& v1,
                                                     const Vector3& v2,
                                                     const Vector3& v3);
  // Calculate a face normal without normalize, no w-information.
  static Vector3 CalculateBasicFaceNormalWithoutNormalize(const Vector3& v1,
                                                          const Vector3& v2,
                                                          const Vector3& v3);

  // Generates a value based on the Gaussian (normal) distribution function
  // with the given offset and scale parameters.
  static Real GaussianDistribution(Real x,
                                   Real offset = 0.0f,
                                   Real scale = 1.0f);

  // Clamp a value within an inclusive range.
  template <typename T>
  static T Clamp(T val, T minval, T maxval) {
    assert(minval <= maxval && "Invalid clamp range");
    return std::max(std::min(val, maxval), minval);
  }

  static Matrix4 MakeViewMatrix(const Vector3& position,
                                const Quaternion& orientation,
                                const Matrix4* reflectMatrix = 0);

  // Get a bounding radius value from a bounding box.
  static Real BoundingRadiusFromAABB(const AxisAlignedBox& aabb);



  static const Real POS_INFINITY;
  static const Real NEG_INFINITY;
  static const Real PI;
  static const Real TWO_PI;
  static const Real HALF_PI;

 protected:
  // Private function to build trig tables.
  void BuildTrigTables();

  static Real SinTable(Real value);
  static Real TanTable(Real value);
  // angle units used by the api
  static AngleUnit angle_unit_;

  /// Size of the trig tables as determined by constructor.
  static int trig_table_size_;

  /// Radian -> index factor value ( mTrigTableSize / 2 * PI )
  static Real trig_table_factor_;
  static Real* sin_table_;
  static Real* tan_table_;

  static const Real deg_to_rad_;
  static const Real rad_to_deg_;
};

}  // namespace math
}  // namespace mist

#endif  // MIST_MATH_MATH_UTIL_H_
