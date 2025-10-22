/*
 * Copyright (C) 2020-2025 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/math/Matrix3.h>
#include <core/math/Vector3.h>

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>

#include <type_traits>

namespace walberla {

template <typename T, typename U = T> inline U es2walberla(T const &v) {
  static_assert(std::is_arithmetic_v<T> and std::is_arithmetic_v<U>,
                "a partial specialization of es2walberla is missing for the "
                "type you are trying to convert");
  return numeric_cast<U>(v);
}
template <> inline Vector3<float> es2walberla(Utils::Vector3d const &v) {
  return Vector3<float>{numeric_cast<float>(v[0]), numeric_cast<float>(v[1]),
                        numeric_cast<float>(v[2])};
}
template <> inline Vector3<double> es2walberla(Utils::Vector3d const &v) {
  return Vector3<double>{v[0], v[1], v[2]};
}

template <typename T> auto to_vector3d(Vector3<T> const &v) noexcept {
  return Utils::Vector3d{Utils::detail::carray_alias<double, 3u>{
      double_c(v[0]), double_c(v[1]), double_c(v[2])}};
}

inline Utils::Vector3i to_vector3i(Vector3<int> const &v) noexcept {
  return {v[0], v[1], v[2]};
}

inline Utils::Vector3i to_vector3i(Cell const &v) noexcept {
  return {v.x(), v.y(), v.z()};
}

template <typename T> auto to_vector3(Utils::Vector3d const &v) noexcept {
  return Vector3<T>{numeric_cast<T>(v[0]), numeric_cast<T>(v[1]),
                    numeric_cast<T>(v[2])};
}

template <typename T> auto to_vector9d(Matrix3<T> const &m) noexcept {
  return Utils::VectorXd<9>{Utils::detail::carray_alias<double, 9u>{
      double_c(m[0]), double_c(m[1]), double_c(m[2]), double_c(m[3]),
      double_c(m[4]), double_c(m[5]), double_c(m[6]), double_c(m[7]),
      double_c(m[8])}};
}

template <typename Function>
void interpolate_bspline_at_pos(Utils::Vector3d const &pos, Function const &f) {
  Utils::Interpolation::bspline_3d<2>(
      pos, f, Utils::Vector3d::broadcast(1.), // grid spacing
      Utils::Vector3d::broadcast(.5));        // offset
}

} // namespace walberla
