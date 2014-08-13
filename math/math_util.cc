// Copyright (c) 2014, Yan, Hongzhang. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include "math/math_util.h"

namespace mist {
namespace math {

const Real Math::POS_INFINITY = std::numeric_limits<Real>::infinity();
const Real Math::NEG_INFINITY = -std::numeric_limits<Real>::infinity();
const Real Math::PI = Real(4.0 * atan(1.0));
const Real Math::TWO_PI = Real(2.0 * PI);
const Real Math::HALF_PI = Real(0.5 * PI);
const Real Math::deg_to_rad_ = PI / Real(180.0);
const Real Math::rad_to_deg_ = Real(180.0) / PI;
const Real Math::LOG2 = log(Real(2.0));

int Math::trig_table_size_;
Math::AngleUnit Math::angle_unit_;

Real Math::trig_table_factor_;
Real* Math::sin_table_ = NULL;
Real* Math::tan_table_ = NULL;

Math::Math(unsigned int trig_table_size) {
  angle_unit_ = AU_DEGREE;
  trig_table_size_ = trig_table_size;
  trig_table_factor_ = trig_table_size_ / Math::TWO_PI;

  sin_table_ = OGRE_ALLOC_T(Real, trig_table_size_, MEMCATEGORY_GENERAL);
  tan_table_ = OGRE_ALLOC_T(Real, trig_table_size_, MEMCATEGORY_GENERAL);

  buildTrigTables();
}

Math::~Math() {
  OGRE_FREE(sin_table_, MEMCATEGORY_GENERAL);
  OGRE_FREE(tan_table_, MEMCATEGORY_GENERAL);
}

void Math::BuildTrigTables(void) {
  // Build trig lookup tables
  // Could get away with building only PI sized Sin table but simpler this
  // way. Who cares, it'll ony use an extra 8k of memory anyway and I like
  // simplicity.
  Real angle;
  for (int i = 0; i < trig_table_size_; ++i) {
    angle = Math::TWO_PI * i / trig_table_size_;
    sin_table_[i] = sin(angle);
    tan_table_[i] = tan(angle);
  }
}

Real Math::SinTable(Real value) {
  // Convert range to index values, wrap if required
  int idx;
  if (value >= 0) {
      idx = static_cast<int>(value * trig_table_factor_) % trig_table_size_;
  } else {
      idx = trig_table_size_ -
            (static_cast<int>(-value * trig_table_factor_) % trig_table_size_) -
            1;
  }

  return sin_table_[idx];
}

Real Math::TanTable(Real value) {
  // Convert range to index values, wrap if required
  int idx = static_cast<int>(value *= trig_table_factor_) % trig_table_size_;
  return tan_table_[idx];
}


int Math::ISign(int value) {
  return ( value > 0 ? +1 : ( value < 0 ? -1 : 0 ) );
}

Real Math::Sign(Real value) {
  if (value > 0.0)
    return 1.0;

  if (value < 0.0)
    return -1.0;
  return 0.0;
}

Radian Math::ACos(Real value) {
  if (-1.0 < value) {
    if (value < 1.0)
      return Radian(acos(value));
    else
      return Radian(0.0);
  } else {
    return Radian(PI);
  }
}

Radian Math::ASin(Real value) {
  if (-1.0 < value) {
    if (value < 1.0)
      return Radian(asin(value));
    else
      return Radian(HALF_PI);
  } else {
    return Radian(-HALF_PI);
  }
}

Real Math::InvSqrt(Real value) {
  return Real(asm_rsq(value));
}

Real Math::UnitRandom() {
    return asm_rand() / asm_rand_max();
}


Real Math::RangeRandom(Real low, Real high) {
  return (high - low)*UnitRandom() + low;
}


Real Math::SymmetricRandom() {
  return 2.0f * UnitRandom() - 1.0f;
}


void Math::SetAngleUnit(Math::AngleUnit unit) {
  angle_unit_ = unit;
}

Math::AngleUnit Math::getAngleUnit(void) {
  return angle_unit_;
}

Real Math::AngleUnitsToRadians(Real angleunits) {
  if (angle_unit_ == AU_DEGREE)
    return angleunits * deg_to_rad_;
  else
    return angleunits;
}

Real Math::RadiansToAngleUnits(Real radians) {
  if (angle_unit_ == AU_DEGREE)
    return radians * rad_to_deg_;
  else
    return radians;
}


Real Math::AngleUnitsToDegrees(Real angleunits) {
  if (angle_unit_ == AU_RADIAN)
    return angleunits * rad_to_deg_;
  else
    return angleunits;
}


Real Math::DegreesToAngleUnits(Real degrees) {
  if (angle_unit_ == AU_RADIAN)
    return degrees * deg_to_rad_;
  else
    return degrees;
}


bool Math::PointInTri2D(const Vector2& p,
                        const Vector2& a,
                        const Vector2& b,
                        const Vector2& c) {
  // Winding must be consistent from all edges for point to be inside
  Vector2 v1, v2;
  Real dot[3];
  bool zero_dot[3];

  v1 = b - a;
  v2 = p - a;
  // Note we don't care about normalisation here since sign is all we need
  // It means we don't have to worry about magnitude of cross products either
  dot[0] = v1.CrossProduct(v2);
  zero_dot[0] = Math::RealEqual(dot[0], 0.0f, 1e-3);

  v1 = c - b;
  v2 = p - b;

  dot[1] = v1.CrossProduct(v2);
  zero_dot[1] = Math::RealEqual(dot[1], 0.0f, 1e-3);

  // Compare signs (ignore colinear / coincident points)
  if (!zero_dot[0] && !zero_dot[1] &&
      Math::Sign(dot[0]) != Math::Sign(dot[1])) {
    return false;
  }

  v1 = a - c;
  v2 = p - c;

  dot[2] = v1.CrossProduct(v2);
  zero_dot[2] = Math::RealEqual(dot[2], 0.0f, 1e-3);
  // Compare signs (ignore colinear / coincident points)
  if (!zero_dot[0] && !zero_dot[2] && Math::Sign(dot[0]) != Math::Sign(dot[2]))
    return false;
  if (!zero_dot[1] && !zero_dot[2] && Math::Sign(dot[1]) != Math::Sign(dot[2]))
    return false;

  return true;
}

bool Math::PointInTri3D(const Vector3& p,
                        const Vector3& a,
                        const Vector3& b,
                        const Vector3& c,
                        const Vector3& normal) {
  // Winding must be consistent from all edges for point to be inside
  Vector3 v1, v2;
  Real dot[3];
  bool zero_dot[3];

  v1 = b - a;
  v2 = p - a;

  // Note we don't care about normalisation here since sign is all we need
  // It means we don't have to worry about magnitude of cross products either
  dot[0] = v1.CrossProduct(v2).DotProduct(normal);
  zero_dot[0] = Math::RealEqual(dot[0], 0.0f, 1e-3);

  v1 = c - b;
  v2 = p - b;

  dot[1] = v1.CrossProduct(v2).DotProduct(normal);
  zero_dot[1] = Math::RealEqual(dot[1], 0.0f, 1e-3);

  // Compare signs (ignore colinear / coincident points)
  if (!zero_dot[0] && !zero_dot[1] &&
      Math::Sign(dot[0]) != Math::Sign(dot[1])) {
          return false;
  }

  v1 = a - c;
  v2 = p - c;

  dot[2] = v1.CrossProduct(v2).DotProduct(normal);
  zero_dot[2] = Math::RealEqual(dot[2], 0.0f, 1e-3);
  // Compare signs (ignore colinear / coincident points)
  if (!zero_dot[0] && !zero_dot[2] && Math::Sign(dot[0]) != Math::Sign(dot[2]))
    return false;
  if (!zero_dot[1] && !zero_dot[2] && Math::Sign(dot[1]) != Math::Sign(dot[2]))
    return false;

  return true;
}

bool Math::RealEqual(Real a, Real b, Real tolerance) {
  if (fabs(b-a) <= tolerance)
    return true;
  else
    return false;
}


std::pair<bool, Real> Math::Intersects(const Ray& ray, const Plane& plane) {
  Real denom = plane.normal.DotProduct(ray.getDirection());
  if (Math::Abs(denom) < std::numeric_limits<Real>::epsilon()) {
      // Parallel
      return std::pair<bool, Real>(false, 0);
  } else {
      Real nom = plane.normal.DotProduct(ray.getOrigin()) + plane.d;
      Real t = -(nom/denom);
      return std::pair<bool, Real>(t >= 0, t);
  }
}

std::pair<bool, Real> Math::Intersects(const Ray& ray,
                                       const vector<Plane>::type& planes,
                                       bool normal_is_outside) {
  list<Plane>::type planesList;
  vector<Plane>::type::const_iterator i = planes.begin();
  for (; i != planes.end(); ++i) {
    planesList.push_back(*i);
  }

  return Intersects(ray, planesList, normal_is_outside);
}

std::pair<bool, Real> Math::Intersects(const Ray& ray,
                                       const list<Plane>::type& planes,
                                       bool normal_is_outside) {
  list<Plane>::type::const_iterator planeit, planeitend;
  planeitend = planes.end();
  bool allInside = true;
  std::pair<bool, Real> ret;
  std::pair<bool, Real> end;
  ret.first = false;
  ret.second = 0.0f;
  end.first = false;
  end.second = 0;

  // derive side
  // NB we don't pass directly since that would require Plane::Side in
  // interface, which results in recursive includes since Math is so fundamental
  Plane::Side outside =
      normal_is_outside ? Plane::POSITIVE_SIDE : Plane::NEGATIVE_SIDE;

  for (planeit = planes.begin(); planeit != planeitend; ++planeit) {
    const Plane& plane = *planeit;
    // is origin outside?
    if (plane.getSide(ray.getOrigin()) == outside) {
      allInside = false;
      // Test single plane
      std::pair<bool, Real> planeRes = ray.Intersects(plane);
      if (planeRes.first) {
        // Ok, we intersected
        ret.first = true;
        // Use the most distant result since convex volume
        ret.second = std::max(ret.second, planeRes.second);
      } else {
        ret.first = false;
        ret.second = 0.0f;
        return ret;
      }
    } else {
      std::pair<bool, Real> planeRes = ray.Intersects(plane);
      if (planeRes.first) {
        if (!end.first) {
          end.first = true;
          end.second = planeRes.second;
        } else {
          end.second = std::min(planeRes.second, end.second);
        }
      }
    }
  }

  if (allInside) {
    // Intersecting at 0 distance since inside the volume!
    ret.first = true;
    ret.second = 0.0f;
    return ret;
  }

  if (end.first) {
    if (end.second < ret.second) {
      ret.first = false;
      return ret;
    }
  }
  return ret;
}

std::pair<bool, Real> Math::Intersects(const Ray& ray,
                                       const Sphere& sphere,
                                       bool discardInside) {
  const Vector3& raydir = ray.getDirection();
  // Adjust ray origin relative to sphere center
  const Vector3& rayorig = ray.getOrigin() - sphere.getCenter();
  Real radius = sphere.getRadius();

  // Check origin inside first
  if (rayorig.squaredLength() <= radius*radius && discardInside) {
      return std::pair<bool, Real>(true, 0);
  }

  // Mmm, quadratics
  // Build coeffs which can be used with std quadratic solver
  // ie t = (-b +/- sqrt(b*b + 4ac)) / 2a
  Real a = raydir.DotProduct(raydir);
  Real b = 2 * rayorig.DotProduct(raydir);
  Real c = rayorig.DotProduct(rayorig) - radius*radius;

  // Calc determinant
  Real d = (b*b) - (4 * a * c);
  if (d < 0) {
    // No intersection
    return std::pair<bool, Real>(false, 0);
  } else {
    // BTW, if d=0 there is one intersection, if d > 0 there are 2
    // But we only want the closest one, so that's ok, just use the
    // '-' version of the solver
    Real t = (-b - Math::Sqrt(d)) / (2 * a);
    if (t < 0)
      t = (-b + Math::Sqrt(d)) / (2 * a);
    return std::pair<bool, Real>(true, t);
  }
}

std::pair<bool, Real> Math::Intersects(const Ray& ray,
                                       const AxisAlignedBox& box) {
  if (box.isNull()) return std::pair<bool, Real>(false, 0);
  if (box.isInfinite()) return std::pair<bool, Real>(true, 0);

  Real lowt = 0.0f;
  Real t;
  bool hit = false;
  Vector3 hitpoint;
  const Vector3& min = box.getMinimum();
  const Vector3& max = box.getMaximum();
  const Vector3& rayorig = ray.getOrigin();
  const Vector3& raydir = ray.getDirection();

  // Check origin inside first
  if (rayorig > min && rayorig < max) {
    return std::pair<bool, Real>(true, 0);
  }

  // Check each face in turn, only check closest 3
  // Min x
  if (rayorig.x <= min.x && raydir.x > 0) {
    t = (min.x - rayorig.x) / raydir.x;
    if (t >= 0) {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.y >= min.y && hitpoint.y <= max.y &&
          hitpoint.z >= min.z && hitpoint.z <= max.z &&
          (!hit || t < lowt)) {
        hit = true;
        lowt = t;
      }
    }
  }
  // Max x
  if (rayorig.x >= max.x && raydir.x < 0) {
    t = (max.x - rayorig.x) / raydir.x;
    if (t >= 0) {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.y >= min.y && hitpoint.y <= max.y &&
          hitpoint.z >= min.z && hitpoint.z <= max.z &&
          (!hit || t < lowt)) {
        hit = true;
        lowt = t;
      }
    }
  }
  // Min y
  if (rayorig.y <= min.y && raydir.y > 0) {
    t = (min.y - rayorig.y) / raydir.y;
    if (t >= 0) {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
          hitpoint.z >= min.z && hitpoint.z <= max.z &&
          (!hit || t < lowt)) {
        hit = true;
        lowt = t;
      }
    }
  }
  // Max y
  if (rayorig.y >= max.y && raydir.y < 0) {
    t = (max.y - rayorig.y) / raydir.y;
    if (t >= 0) {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
          hitpoint.z >= min.z && hitpoint.z <= max.z &&
          (!hit || t < lowt)) {
        hit = true;
        lowt = t;
      }
    }
  }
  // Min z
  if (rayorig.z <= min.z && raydir.z > 0) {
    t = (min.z - rayorig.z) / raydir.z;
    if (t >= 0) {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
          hitpoint.y >= min.y && hitpoint.y <= max.y &&
          (!hit || t < lowt)) {
        hit = true;
        lowt = t;
      }
    }
  }
  // Max z
  if (rayorig.z >= max.z && raydir.z < 0) {
    t = (max.z - rayorig.z) / raydir.z;
    if (t >= 0) {
      // Substitute t back into ray and check bounds and dist
      hitpoint = rayorig + raydir * t;
      if (hitpoint.x >= min.x && hitpoint.x <= max.x &&
          hitpoint.y >= min.y && hitpoint.y <= max.y &&
          (!hit || t < lowt)) {
        hit = true;
        lowt = t;
      }
    }
  }

  return std::pair<bool, Real>(hit, lowt);
}

bool Math::Intersects(const Ray& ray,
                      const AxisAlignedBox& box,
                      Real* d1,
                      Real* d2) {
  if (box.isNull())
    return false;

  if (box.isInfinite()) {
    if (d1) *d1 = 0;
    if (d2) *d2 = Math::POS_INFINITY;
    return true;
  }

  const Vector3& min = box.getMinimum();
  const Vector3& max = box.getMaximum();
  const Vector3& rayorig = ray.getOrigin();
  const Vector3& raydir = ray.getDirection();

  Vector3 absDir;
  absDir[0] = Math::Abs(raydir[0]);
  absDir[1] = Math::Abs(raydir[1]);
  absDir[2] = Math::Abs(raydir[2]);

  // Sort the axis, ensure check minimise floating error axis first
  int imax = 0, imid = 1, imin = 2;
  if (absDir[0] < absDir[2]) {
    imax = 2;
    imin = 0;
  }

  if (absDir[1] < absDir[imin]) {
    imid = imin;
    imin = 1;
  } else if (absDir[1] > absDir[imax]) {
    imid = imax;
    imax = 1;
  }

  Real start = 0, end = Math::POS_INFINITY;

#define _CALC_AXIS(i)                                     \
  do {                                                    \
      Real denom = 1 / raydir[i];                         \
      Real newstart = (min[i] - rayorig[i]) * denom;      \
      Real newend = (max[i] - rayorig[i]) * denom;        \
      if (newstart > newend) std::swap(newstart, newend); \
      if (newstart > end || newend < start) return false; \
      if (newstart > start) start = newstart;             \
      if (newend < end) end = newend;                     \
  } while (0)

    // Check each axis in turn

  _CALC_AXIS(imax);

  if (absDir[imid] < std::numeric_limits<Real>::epsilon()) {
      // Parallel with middle and minimise axis, check bounds only
      if (rayorig[imid] < min[imid] || rayorig[imid] > max[imid] ||
          rayorig[imin] < min[imin] || rayorig[imin] > max[imin])
          return false;
  } else {
    _CALC_AXIS(imid);

    if (absDir[imin] < std::numeric_limits<Real>::epsilon()) {
      // Parallel with minimise axis, check bounds only
      if (rayorig[imin] < min[imin] || rayorig[imin] > max[imin])
        return false;
    } else {
      _CALC_AXIS(imin);
    }
  }
#undef _CALC_AXIS

  if (d1) *d1 = start;
  if (d2) *d2 = end;

  return true;
}

std::pair<bool, Real> Math::Intersects(const Ray& ray,
                                       const Vector3& a,
                                       const Vector3& b,
                                       const Vector3& c,
                                       const Vector3& normal,
                                       bool positive_side,
                                       bool negative_side) {
  // Calculate intersection with plane.
  Real t;
  {
    Real denom = normal.DotProduct(ray.getDirection());

    // Check intersect side
    if (denom > + std::numeric_limits<Real>::epsilon()) {
      if (!negative_side)
        return std::pair<bool, Real>(false, 0);
    } else if (denom < - std::numeric_limits<Real>::epsilon()) {
      if (!positive_side)
        return std::pair<bool, Real>(false, 0);
    } else {
      // Parallel or triangle area is close to zero when
      // the plane normal not Normalised.
      return std::pair<bool, Real>(false, 0);
    }

    t = normal.DotProduct(a - ray.getOrigin()) / denom;

    if (t < 0) {
      // Intersection is behind origin
      return std::pair<bool, Real>(false, 0);
    }
  }

  // Calculate the largest area projection plane in X, Y or Z.
  size_t i0, i1;
  {
    Real n0 = Math::Abs(normal[0]);
    Real n1 = Math::Abs(normal[1]);
    Real n2 = Math::Abs(normal[2]);

    i0 = 1; i1 = 2;
    if (n1 > n2) {
      if (n1 > n0)
        i0 = 0;
    } else {
      if (n2 > n0)
        i1 = 0;
    }
  }

  // Check the intersection point is inside the triangle.
  {
    Real u1 = b[i0] - a[i0];
    Real v1 = b[i1] - a[i1];
    Real u2 = c[i0] - a[i0];
    Real v2 = c[i1] - a[i1];
    Real u0 = t * ray.getDirection()[i0] + ray.getOrigin()[i0] - a[i0];
    Real v0 = t * ray.getDirection()[i1] + ray.getOrigin()[i1] - a[i1];

    Real alpha = u0 * v2 - u2 * v0;
    Real beta  = u1 * v0 - u0 * v1;
    Real area  = u1 * v2 - u2 * v1;

    // epsilon to avoid float precision error
    const Real EPSILON = 1e-6f;

    Real tolerance = - EPSILON * area;

    if (area > 0) {
      if (alpha < tolerance || beta < tolerance || alpha+beta > area-tolerance)
        return std::pair<bool, Real>(false, 0);
    } else {
      if (alpha > tolerance || beta > tolerance || alpha+beta < area-tolerance)
        return std::pair<bool, Real>(false, 0);
    }
  }

  return std::pair<bool, Real>(true, t);
}

std::pair<bool, Real> Math::Intersects(const Ray& ray,
                                       const Vector3& a,
                                       const Vector3& b,
                                       const Vector3& c,
                                       bool positive_side,
                                       bool negative_side) {
  Vector3 normal = CalculateBasicFaceNormalWithoutNormalize(a, b, c);
  return Intersects(ray, a, b, c, normal, positive_side, negative_side);
}

bool Math::Intersects(const Sphere& sphere, const AxisAlignedBox& box) {
  if (box.isNull())
    return false;
  if (box.isInfinite())
    return true;

  // Use splitting planes
  const Vector3& center = sphere.getCenter();
  Real radius = sphere.getRadius();
  const Vector3& min = box.getMinimum();
  const Vector3& max = box.getMaximum();

  // Arvo's algorithm
  Real s, d = 0;
  for (int i = 0; i < 3; ++i) {
    if (center.ptr()[i] < min.ptr()[i]) {
      s = center.ptr()[i] - min.ptr()[i];
      d += s * s;
    } else if (center.ptr()[i] > max.ptr()[i]) {
      s = center.ptr()[i] - max.ptr()[i];
      d += s * s;
    }
  }

  return d <= radius * radius;
}

bool Math::Intersects(const Plane& plane, const AxisAlignedBox& box) {
  return (plane.getSide(box) == Plane::BOTH_SIDE);
}

bool Math::Intersects(const Sphere& sphere, const Plane& plane) {
  return
      (Math::Abs(plane.getDistance(sphere.getCenter())) <= sphere.getRadius() );
}

Vector3 Math::CalculateTangentSpaceVector(const Vector3& position1,
                                          const Vector3& position2,
                                          const Vector3& position3,
                                          Real u1, Real v1,
                                          Real u2, Real v2,
                                          Real u3, Real v3) {
  // side0 is the vector along one side of the triangle of vertices passed in,
  // and side1 is the vector along another side. Taking the cross product of
  // these returns the normal.
  Vector3 side0 = position1 - position2;
  Vector3 side1 = position3 - position1;
  // Calculate face normal
  Vector3 normal = side1.CrossProduct(side0);
  normal.Normalise();
  // Now we use a formula to Calculate the tangent.
  Real deltaV0 = v1 - v2;
  Real deltaV1 = v3 - v1;
  Vector3 tangent = deltaV1 * side0 - deltaV0 * side1;
  tangent.Normalise();
  // Calculate binormal
  Real deltaU0 = u1 - u2;
  Real deltaU1 = u3 - u1;
  Vector3 binormal = deltaU1 * side0 - deltaU0 * side1;
  binormal.Normalise();
  // Now, we take the cross product of the tangents to get a vector which
  // should point in the same direction as our normal Calculated above.
  // If it points in the opposite direction (the dot product between the normals
  // is less than zero), then we need to reverse the s and t tangents.
  // This is because the triangle has been mirrored when going from tangent
  // space to object space. reverse tangents if necessary.
  Vector3 tangentCross = tangent.CrossProduct(binormal);
  if (tangentCross.DotProduct(normal) < 0.0f) {
    tangent = -tangent;
    binormal = -binormal;
  }

  return tangent;
}

Matrix4 Math::BuildReflectionMatrix(const Plane& p) {
  return Matrix4(
      -2 * p.normal.x * p.normal.x + 1,
      -2 * p.normal.x * p.normal.y,
      -2 * p.normal.x * p.normal.z,
      -2 * p.normal.x * p.d,

      -2 * p.normal.y * p.normal.x,
      -2 * p.normal.y * p.normal.y + 1,
      -2 * p.normal.y * p.normal.z,
      -2 * p.normal.y * p.d,

      -2 * p.normal.z * p.normal.x,
      -2 * p.normal.z * p.normal.y,
      -2 * p.normal.z * p.normal.z + 1,
      -2 * p.normal.z * p.d,

      0, 0, 0, 1);
}

Vector4 Math::CalculateFaceNormal(const Vector3& v1,
                                  const Vector3& v2,
                                  const Vector3& v3) {
  Vector3 normal = CalculateBasicFaceNormal(v1, v2, v3);
  // Now set up the w (distance of tri from origin
  return Vector4(normal.x, normal.y, normal.z, -(normal.DotProduct(v1)));
}

Vector3 Math::CalculateBasicFaceNormal(const Vector3& v1,
                                       const Vector3& v2,
                                       const Vector3& v3) {
  Vector3 normal = (v2 - v1).CrossProduct(v3 - v1);
  normal.Normalise();
  return normal;
}

Vector4 Math::CalculateFaceNormalWithoutNormalize(const Vector3& v1,
                                                  const Vector3& v2,
                                                  const Vector3& v3) {
  Vector3 normal = CalculateBasicFaceNormalWithoutNormalize(v1, v2, v3);
  // Now set up the w (distance of tri from origin)
  return Vector4(normal.x, normal.y, normal.z, -(normal.DotProduct(v1)));
}

Vector3 Math::CalculateBasicFaceNormalWithoutNormalize(const Vector3& v1,
                                                       const Vector3& v2,
                                                       const Vector3& v3) {
  Vector3 normal = (v2 - v1).CrossProduct(v3 - v1);
  return normal;
}

Real Math::GaussianDistribution(Real x, Real offset, Real scale) {
  Real nom = Math::Exp(-Math::Sqr(x - offset) / (2 * Math::Sqr(scale)));
  Real denom = scale * Math::Sqrt(2 * Math::PI);

  return nom / denom;
}

Matrix4 Math::MakeViewMatrix(const Vector3& position,
                             const Quaternion& orientation,
                             const Matrix4* reflect_matrix) {
  Matrix4 view_matrix;

  // View matrix is:
  //
  //  [ Lx  Uy  Dz  Tx  ]
  //  [ Lx  Uy  Dz  Ty  ]
  //  [ Lx  Uy  Dz  Tz  ]
  //  [ 0   0   0   1   ]
  //
  // Where T = -(Transposed(Rot) * Pos)

  // This is most efficiently done using 3x3 Matrices
  Matrix3 rot;
  orientation.ToRotationMatrix(rot);

  // Make the translation relative to new axes
  Matrix3 rotT = rot.Transpose();
  Vector3 trans = -rotT * position;

  // Make final matrix
  view_matrix = Matrix4::IDENTITY;
  // fills upper 3x3
  view_matrix = rotT;
  view_matrix[0][3] = trans.x;
  view_matrix[1][3] = trans.y;
  view_matrix[2][3] = trans.z;

  // Deal with reflections
  if (reflect_matrix) {
    view_matrix = view_matrix * (*reflect_matrix);
  }

  return view_matrix;
}

Real Math::BoundingRadiusFromAABB(const AxisAlignedBox& aabb) {
  Vector3 max = aabb.GetMaximum();
  Vector3 min = aabb.GetMinimum();

  Vector3 magnitude = max;
  magnitude.MakeCeil(-max);
  magnitude.MakeCeil(min);
  magnitude.MakeCeil(-min);

  return magnitude.Length();
}

}  // namespace math
}  // namespace mist

