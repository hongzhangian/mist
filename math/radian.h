// Copyright (c) 2014, Yan, Hongzhang. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#ifndef MIST_MATH_RADIAN_H_
#define MIST_MATH_RADIAN_H_

#include "math/math_util.h"

namespace mist {
namespace math {

class Degree;

class Radian {
 public:
  inline explicit Radian(Real rad = 0)
      : rad_(rad) {
  }

  // these functions could not be defined within the class definition of class
  // Radian because they required class Degree to be defined
  inline explicit Radian(const Degree& degree)
      : rad_(degree.ValueAsRadian()) {
  }

  inline Radian ValueAsRadian() const { return rad_; }
  inline Degree ValueAsDegrees() const  {
    return Math::RadiansToDegrees(rad_);
  }
  inline Real ValueAsAngleUnits() const {
    return Math::RadiansToAngleUnits (rad_);
  }

  inline Radian& operator=(const Real& real) {
    rad_ = real;
    return *this;
  }
  inline Radian& operator=(const Radian& radian) {
    rad_ = radian.rad_;
    return *this;
  }
  inline Radian& operator=(const Degree& degree) {
    rad_ = degree.AsRadians();
    return *this;
  }

  inline const Radian& operator+() const { return *this; }
  inline Radian operator+(const Radian& radian) const {
    return Radian (rad_ + radian.rad_);
  }
  inline Radian Radian::operator+(const Degree& degree) const {
    return Radian ( rad_ + degree.AsRadians() );
  }

  inline Radian& operator+=(const Radian& radian) {
    rad_ += radian.rad_;
    return *this;
  }
  inline Radian& Radian::operator+=(const Degree& degree) {
    rad_ += degree.AsRadians();
    return *this;
  }

  inline Radian operator-() const { return Radian(-rad_); }
  inline Radian operator-(const Radian& radian) const {
    return Radian ( rad_ - radian.rad_ );
  }
  inline Radian operator-(const Degree& degree) const {
    return Radian ( rad_ - degree.AsRadians() );
  }

  inline Radian& operator-= (const Radian& radian) {
    rad_ -= radian.rad_;
    return *this;
  }
  inline Radian& Radian::operator-=(const Degree& degree) {
    rad_ -= degree.AsRadians();
    return *this;
  }

  inline Radian operator*(Real real) const {
    return Radian ( rad_ * real);
  }
  inline Radian operator*(const Radian& radian) const {
    return Radian ( rad_ * radian.rad_ );
  }
  inline friend Radian operator*(Real a, const Radian& b) {
    return Radian (a * b.ValueAsRadian());
  }
  inline Radian& operator*=( Real real ) { rad_ *= real; return *this; }

  inline Radian operator/(Real real) const {
    return Radian ( rad_ / real);
  }
  inline Radian operator/(const Radian& radian) const {
    return Radian ( rad_ / radian.rad_ );
  }
  inline friend Radian operator/(Real a, const Radian& b) {
    return Radian (a / b.ValueAsRadian());
  }
  inline Radian& operator/=(Real real) {
    rad_ /= real;
    return *this;
  }

  bool operator< (const Radian& radian) const { return rad_ <  radian.rad_; }
  bool operator<=(const Radian& radian) const { return rad_ <= radian.rad_; }
  bool operator==(const Radian& radian) const { return rad_ == radian.rad_; }
  bool operator!=(const Radian& radian) const { return rad_ != radian.rad_; }
  bool operator>=(const Radian& radian) const { return rad_ >= radian.rad_; }
  bool operator> (const Radian& radian) const { return rad_ >  radian.rad_; }

 private:
  Real rad_;
};

}  // namespace math
}  // namespace mist

#endif  // MIST_MATH_RADIAN_H_
