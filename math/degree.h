// Copyright (c) 2014, Yan, Hongzhang. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#ifndef MIST_MATH_DEGREE_H_
#define MIST_MATH_DEGREE_H_

#include "math/math_util.h"

namespace mist {
namespace math {

class Degree {
 public:
  explicit Degree(Real deg = 0)
      : deg_(deg) {
  }

  explicit Degree(const Radian& radian)
      : deg_(radian.ValueAsDegrees()) {
  }


  inline Real ValueAsDegrees() const { return deg_; }
  inline Real ValueAsRadians() const {
    return Math::DegreesToRadians (deg_);
  }
  inline Real ValueAsAngleUnits() const {
    return Math::DegreesToAngleUnits (deg_);
  }


  inline Degree& operator=(const Real& real) {
    deg_ = real;
    return *this;
  }
  inline Degree& operator=(const Degree& degree) {
    deg_ = degree.deg_;
    return *this;
  }
  inline Degree& operator=(const Radian& radian) {
    deg_ = radian.ValueAsDegrees();
    return *this;
  }

  inline const Degree& operator+() const { return *this; }
  inline Degree operator+(const Degree& degree) const {
    return Degree (deg_ + degree.deg_);
  }
  inline Degree operator+(const Radian& radian) const {
    return Degree (deg_ + radian.ValueAsDegrees());
  }
  inline Degree& operator+=(const Degree& degree) {
    deg_ += degree.deg_;
    return *this;
  }
  inline Degree& operator+=(const Radian& radian) {
    deg_ += radian.ValueAsDegrees();
    return *this;
  }

  inline Degree operator-() const { return Degree(-deg_); }
  inline Degree operator-(const Degree& degree) const {
    return Degree (deg_ - degree.deg_);
  }
  inline Degree operator-(const Radian& radian) const {
    return Degree (deg_ - radian.ValueAsDegrees());
  }
  inline Degree& operator-=(const Degree& degree) {
    deg_ -= degree.deg_; return *this;
  }
  inline Degree& operator-=(const Radian& radian) {
    deg_ -= radian.ValueAsDegrees(); return *this;
  }

  inline Degree operator*(Real real) const {
    return Degree ( deg_ * real );
  }
  inline Degree operator*(const Degree& real) const {
    return Degree ( deg_ * real.deg_ );
  }
  inline friend Degree operator*(Real real, const Degree& degree) {
    return Degree ( real * degree.ValueAsDegrees() );
  }
  inline Degree& operator*=(Real real ) {
    deg_ *= real;
    return *this;
  }

  inline Degree operator/(Real real) const {
    return Degree ( deg_ / real );
  }
  inline friend Degree operator/(Real real, const Degree& degree) {
    return Degree ( real / degree.ValueAsDegrees() );
  }
  inline Degree& operato/=(Real real) {
    deg_ /= real;
    return *this;
  }

  bool operator< (const Degree& degree) const { return deg_ <  degree.deg_; }
  bool operator<=(const Degree& degree) const { return deg_ <= degree.deg_; }
  bool operator==(const Degree& degree) const { return deg_ == degree.deg_; }
  bool operator!=(const Degree& degree) const { return deg_ != degree.deg_; }
  bool operator>=(const Degree& degree) const { return deg_ >= degree.deg_; }
  bool operator> (const Degree& degree) const { return deg_ >  degree.deg_; }

 private:
  Real deg_;
};

}  // namespace math
}  // namespace mist

#endif  // MIST_MATH_DEGREE_H_
