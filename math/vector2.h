// Copyright (c) 2014, Yan, Hongzhang. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#ifndef MIST_MATH_VECTOR2_H_
#define MIST_MATH_VECTOR2_H_

#include <assert.h>

#include <algorithm>

#include "math/radian.h"

namespace mist {
namespace math {

typedef float Real;

class Vector2 {
 private:
  Real x_, y_;


 public:
  inline Vector2() {}

  inline Vector2(const Real x, const Real y)
      : x_(x),
        y_(y) {
  }

  inline explicit Vector2(const Real scaler)
      : x_(scaler),
        y_(scaler) {
  }

  inline explicit Vector2(const Real af_coordinate[2])
    : x_(af_coordinate[0]),
      y_(af_coordinate[1]) {
  }

  inline explicit Vector2(const int af_coordinate[2]) {
    x_ = (Real)af_coordinate[0];
    y_ = (Real)af_coordinate[1];
  }

  inline explicit Vector2(Real* const r)
    : x_(r[0]),
      y_(r[1]) {
  }

  inline Real x() { return x_; }
  inline Real y() { return y_; }

  // Swap the contents of this vector with another.
  inline void swap(Vector2& vector2) {
    std::swap(x_, vector2.x_);
    std::swap(y_, vector2.y_);
  }

  inline Real operator[](const size_t i) const {
    assert(i < 2);
    return *(&x_ + i);
  }

  inline Real& operator[](const size_t i) {
    assert(i < 2);
    return *(&x_ + i);
  }

  // Pointer accessor for direct copying
  inline Real* ptr() {
    return &x_;
  }
  // Pointer accessor for direct copying
  inline const Real* ptr() const {
    return &x_;
  }


  inline Vector2& operator=(const Vector2& vector2) {
    x_ = vector2.x_;
    y_ = vector2.y_;
    return *this;
  }

  inline Vector2& operator=(const Real scalar) {
    x_ = scalar;
    y_ = scalar;
    return *this;
  }

  inline bool operator==(const Vector2& vector2) const {
    return (x_ == vector2.x_ && y_ == vector2.y_);
  }

  inline bool operator!=(const Vector2& vector2) const {
    return (x_ != vector2.x_ || y_ != vector2.y_);
  }

  // arithmetic operations
  inline Vector2 operator+(const Vector2& vector2) const {
    return Vector2(x_ + vector2.x_, y_ + vector2.y_);
  }

  inline Vector2 operator-(const Vector2& vector2) const {
    return Vector2(x_ - vector2.x_, y_ - vector2.y_);
  }

  inline Vector2 operator*(const Real scalar) const {
    return Vector2(x_ * scalar, y_ * scalar);
  }

  inline Vector2 operator*(const Vector2& r_vector2) const {
    return Vector2(x_ * r_vector2.x_, y_ * r_vector2.y_);
  }

  inline Vector2 operator/(const Real scalar) const {
    assert(scalar != 0.0);
    Real inv = 1.0f / scalar;
    return Vector2(x_ * inv, y_ * inv);
  }

  inline Vector2 operator/(const Vector2& r_vector2) const {
    return Vector2(x_ / r_vector2.x_, y_ / r_vector2.y_);
  }

  inline const Vector2& operator+() const {
    return *this;
  }

  inline Vector2 operator-() const {
    return Vector2(-x_, -y_);
  }

  // overloaded operators to help Vector2
  inline friend Vector2 operator*(const Real scalar, const Vector2& r_vector2) {
    return Vector2(scalar * r_vector2.x_, scalar * r_vector2.y_);
  }

  inline friend Vector2 operator/(const Real scalar, const Vector2& r_vector2) {
    return Vector2(scalar / r_vector2.x_, scalar / r_vector2.y_);
  }

  inline friend Vector2 operator+(const Vector2& vector2, const Real real) {
    return Vector2(vector2.x_ + real, vector2.y_ + real);
  }

  inline friend Vector2 operator+(const Real real, const Vector2& vector2) {
    return Vector2(real + vector2.x_, real + vector2.y_);
  }

  inline friend Vector2 operator-(const Vector2& vector2, const Real real) {
    return Vector2(vector2.x_ - real, vector2.y_ - real);
  }

  inline friend Vector2 operator-(const Real real, const Vector2& vector2) {
    return Vector2(real - vector2.x_, real - vector2.y_);
  }

  // arithmetic updates
  inline Vector2& operator+=(const Vector2& vector2) {
    x_ += vector2.x_;
    y_ += vector2.y_;

    return *this;
  }

  inline Vector2& operator+=(const Real real) {
    x_ += real;
    y_ += real;

    return *this;
  }

  inline Vector2& operator-=(const Vector2& vector2) {
    x_ -= vector2.x_;
    y_ -= vector2.y_;

    return *this;
  }

  inline Vector2& operator-=(const Real real) {
    x_ -= real;
    y_ -= real;

    return *this;
  }

  inline Vector2& operator*=(const Real scalar) {
    x_ *= scalar;
    y_ *= scalar;

    return *this;
  }

  inline Vector2& operator*=(const Vector2& vector2) {
    x_ *= vector2.x_;
    y_ *= vector2.y_;

    return *this;
  }

  inline Vector2& operator/=(const Real scalar) {
    assert(scalar != 0.0);

    Real inv = 1.0f / scalar;

    x_ *= inv;
    y_ *= inv;

    return *this;
  }

  inline Vector2& operator/=(const Vector2& vector2) {
    x_ /= vector2.x_;
    y_ /= vector2.y_;

    return *this;
  }

  // Returns the length (magnitude) of the vector.
  // @warning
  // This operation requires a square root and is expensive in
  // terms of CPU operations. If you don't need to know the exact
  // length (e.g. for just comparing lengths) use squaredLength()
  // instead.
  inline Real Length() const {
    return Sqrt(x_ * x_ + y_ * y_);
  }

  // Returns the square of the length(magnitude) of the vector.
  // @remarks
  // This  method is for efficiency - calculating the actual
  // length of a vector requires a square root, which is expensive
  // in terms of the operations required. This method returns the
  // square of the length of the vector, i.e. the same as the
  // length but before the square root is taken. Use this if you
  // want to find the longest / shortest vector without incurring
  // the square root.
  inline Real SquaredLength() const {
    return x_ * x_ + y_ * y_;
  }

  // Returns the distance to another vector.
  // @warning
  // This operation requires a square root and is expensive in
  // terms of CPU operations. If you don't need to know the exact
  // distance (e.g. for just comparing distances) use squaredDistance()
  // instead.
  inline Real Distance(const Vector2& vector2) const {
    return (*this - vector2).Length();
  }

  // Returns the square of the distance to another vector.
  // @remarks
  // This method is for efficiency - calculating the actual
  // distance to another vector requires a square root, which is
  // expensive in terms of the operations required. This method
  // returns the square of the distance to another vector, i.e.
  // the same as the distance but before the square root is taken.
  // Use this if you want to find the longest / shortest distance
  // without incurring the square root.
  inline Real SquaredDistance(const Vector2& vector2) const {
    return (*this - vector2).SquaredLength();
  }

  // Calculates the dot (scalar) product of this vector with another.
  // @remarks
  // The dot product can be used to calculate the angle between 2
  // vectors. If both are unit vectors, the dot product is the
  // cosine of the angle; otherwise the dot product must be
  // divided by the product of the lengths of both vectors to get
  // the cosine of the angle. This result can further be used to
  // calculate the distance of a point from a plane.
  // @param
  // vector2 Vector with which to calculate the dot product (together
  // with this one).
  // @return
  // A float representing the dot product value.
  inline Real DotProduct(const Vector2& vector2) const {
    return x_ * vector2.x_ + y_ * vector2.y_;
  }

  // Normalises the vector.
  // @remarks
  // This method normalises the vector such that it's
  // length / magnitude is 1. The result is called a unit vector.
  // @note
  // This function will not crash for zero-sized vectors, but there
  // will be no changes made to their components.
  // @return The previous length of the vector.
  inline Real Normalise() {
    Real length = Math::Sqrt(x_ * x_ + y_ * y_);

    // Will also work for zero-sized vectors, but will change nothing
    // We're not using epsilons because we don't need to.
    // Read http://www.ogre3d.org/forums/viewtopic.php?f=4&t=61259
    if (length > Real(0.0f)) {
      Real inv_length = 1.0f / length;
      x_ *= inv_length;
      y_ *= inv_length;
    }

    return length;
  }

  // Returns a vector at a point half way between this and the passed
  // in vector.
  inline Vector2 MidPoint(const Vector2& vector2) const {
    return Vector2((x_ + vector2.x_) * 0.5f, (y_ + vector2.y_) * 0.5f );
  }

  // Returns true if the vector's scalar components are all greater
  // that the ones of the vector it is compared against.
  inline bool operator<(const Vector2& vector2) const {
    return x_ < vector2.x_ && y_ < vector2.y_;
  }

  // Returns true if the vector's scalar components are all smaller
  // that the ones of the vector it is compared against.
  inline bool operator>(const Vector2& vector2) const {
    return x_ > vector2.x_ && y_ > vector2.y_;
  }

  // Sets this vector's components to the minimum of its own and the
  // ones of the passed in vector.
  // @remarks
  // 'Minimum' in this case means the combination of the lowest
  // value of x, y and z from both vectors. Lowest is taken just
  // numerically, not magnitude, so -1 < 0.
  inline void MakeFloor(const Vector2& vector2) {
    if (vector2.x_ < x_)
      x_ = vector2.x_;
    if (vector2.y_ < y_)
      y_ = vector2.y_;
  }

  // Sets this vector's components to the maximum of its own and the
  // ones of the passed in vector.
  // @remarks
  // 'Maximum' in this case means the combination of the highest
  // value of x, y and z from both vectors. Highest is taken just
  // numerically, not magnitude, so 1 > -3.
  inline void MakeCeil(const Vector2& vector2) {
    if (vector2.x_ > x_)
      x_ = vector2.x_;
    if (vector2.y_ > y_)
      y_ = vector2.y_;
  }

  // Generates a vector perpendicular to this vector (eg an 'up' vector).
  // @remarks
  // This method will return a vector which is perpendicular to this
  // vector. There are an infinite number of possibilities but this
  // method will guarantee to generate one of them. If you need more
  // control you should use the Quaternion class.
  inline Vector2 Perpendicular() const {
    return Vector2 (-y_, x_);
  }

  // Calculates the 2 dimensional cross-product of 2 vectors, which results
  // in a single floating point value which is 2 times the area of the triangle.
  inline Real CrossProduct(const Vector2& vector2) const {
    return x_ * vector2.y_ - y_ * vector2.x_;
  }

  // Generates a new random vector which deviates from this vector by a
  // given angle in a random direction.
  // @remarks
  // This method assumes that the random number generator has already
  // been seeded appropriately.
  // @param
  // angle The angle at which to deviate in radians
  // @param
  // up Any vector perpendicular to this one (which could generated
  // by cross-product of this vector and any other non-colinear
  // vector). If you choose not to provide this the function will
  // derive one on it's own, however if you provide one yourself the
  // function will be faster (this allows you to reuse up vectors if
  // you call this method more than once)
  // @return
  // A random vector which deviates from this vector by angle. This
  // vector will not be normalised, normalise it if you wish
  // afterwards.
  inline Vector2 RandomDeviant(Real angle) const {
    angle *=  Math::UnitRandom() * Math::TWO_PI;
    Real cosa = cos(angle);
    Real sina = sin(angle);
    return  Vector2(cosa * x_ - sina * y_, sina * x_ + cosa * y_);
  }

  // Returns true if this vector is zero length.
  inline bool IsZeroLength() const {
    Real sqlen = (x_ * x_) + (y_ * y_);
    return (sqlen < (1e-06 * 1e-06));
  }

  // As normalise, except that this vector is unaffected and the
  // normalised vector is returned as a copy.
  inline Vector2 NormalisedCopy() const {
    Vector2 ret = *this;
    ret.normalise();
    return ret;
  }

  // Calculates a reflection vector to the plane with the given normal .
  // @remarks NB assumes 'this' is pointing AWAY FROM the plane, invert if it is
  // not.
  inline Vector2 Reflect(const Vector2& normal) const {
    return Vector2(*this - (2 * this->dotProduct(normal) * normal));
  }

  // Check whether this vector contains valid values
  inline bool IsNaN() const {
    return Math::IsNaN(x) || Math::IsNaN(y);
  }

  // Gets the angle between 2 vectors.
  // @remarks
  // Vectors do not have to be unit-length but must represent directions.
  inline Radian AngleBetween(const Vector2& other) const {
    Real len_product = Length() * other.Length();
    // Divide by zero check
    if (len_product < 1e-6f)
      len_product = 1e-6f;

    Real f = DotProduct(other) / len_product;

    f = Clamp(f, (Real)-1.0, (Real)1.0);
    return ACos(f);
  }

  // Gets the oriented angle between 2 vectors.
  // @remarks
  // Vectors do not have to be unit-length but must represent directions.
  // The angle is comprised between 0 and 2 PI.
  inline Radian AngleTo(const Vector2& other) const {
    Radian angle = AngleBetween(other);

    if (CrossProduct(other) < 0)
      angle = (Radian)Math::TWO_PI - angle;

    return angle;
  }

  // special points
  static const Vector2 ZERO;
  static const Vector2 UNIT_X;
  static const Vector2 UNIT_Y;
  static const Vector2 NEGATIVE_UNIT_X;
  static const Vector2 NEGATIVE_UNIT_Y;
  static const Vector2 UNIT_SCALE;
};

}  // namespace math
}  // namespace mist

#endif  // MIST_MATH_VECTOR2_H_
