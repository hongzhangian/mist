// Copyright (c) 2014, Yan, Hongzhang. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include "math/vector3.h"

namespace mist {
namespace math {
  const Vector3 Vector3::ZERO( 0, 0, 0 );

  const Vector3 Vector3::UNIT_X( 1, 0, 0 );
  const Vector3 Vector3::UNIT_Y( 0, 1, 0 );
  const Vector3 Vector3::UNIT_Z( 0, 0, 1 );
  const Vector3 Vector3::NEGATIVE_UNIT_X( -1,  0,  0 );
  const Vector3 Vector3::NEGATIVE_UNIT_Y(  0, -1,  0 );
  const Vector3 Vector3::NEGATIVE_UNIT_Z(  0,  0, -1 );
  const Vector3 Vector3::UNIT_SCALE(1, 1, 1);
}  // namespace math
}  // namespace mist
