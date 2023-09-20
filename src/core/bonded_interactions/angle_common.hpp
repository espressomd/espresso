/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002-2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef ANGLE_COMMON_H
#define ANGLE_COMMON_H
/** \file
 *  Common code for functions calculating angle forces.
 */

#include "config/config.hpp"

#include <utils/Vector.hpp>

#include <cmath>
#include <tuple>

namespace detail {
inline double sanitize_cosine(double cosine) {
  if (cosine > TINY_COS_VALUE)
    cosine = TINY_COS_VALUE;
  if (cosine < -TINY_COS_VALUE)
    cosine = -TINY_COS_VALUE;
  return cosine;
}
} // namespace detail

/** Compute the cosine of the angle between three particles.
 *
 *  @param[in]  vec1             Vector from central particle to left particle.
 *  @param[in]  vec2             Vector from central particle to right particle.
 *  @param[in]  sanitize_cosine  Sanitize the cosine of the angle.
 *  @return @f$ \vec{r_{ij}} @f$, @f$ \vec{r_{kj}} @f$,
 *          @f$ \left\|\vec{r_{ij}}\right\|^{-1} @f$,
 *          @f$ \left\|\vec{r_{kj}}\right\|^{-1} @f$,
 *          @f$ \cos(\theta_{ijk}) @f$
 */
inline double calc_cosine(Utils::Vector3d const &vec1,
                          Utils::Vector3d const &vec2,
                          bool sanitize_cosine = false) {
  /* cosine of the angle between vec1 and vec2 */
  auto cos_phi = (vec1 * vec2) / std::sqrt(vec1.norm2() * vec2.norm2());
  if (sanitize_cosine) {
    cos_phi = detail::sanitize_cosine(cos_phi);
  }
  return cos_phi;
}

/** Compute a three-body angle interaction force.
 *
 *  See the details in @ref bondedIA_angle_force. The @f$ K(\theta_{ijk}) @f$
 *  term is provided as a lambda function in @p forceFactor.
 *
 *  @param[in]  vec1             Vector from central particle to left particle.
 *  @param[in]  vec2             Vector from central particle to right particle.
 *  @param[in]  forceFactor      Angle force term.
 *  @param[in]  sanitize_cosine  Sanitize the cosine of the angle.
 *  @tparam     ForceFactor      Function evaluating the angle force term
 *                               for a given angle.
 *  @return Forces on the second, first and third particles, in that order.
 */
template <typename ForceFactor>
std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
angle_generic_force(Utils::Vector3d const &vec1, Utils::Vector3d const &vec2,
                    ForceFactor forceFactor, bool sanitize_cosine) {
  auto const d1 = vec1.norm();
  auto const d2 = vec2.norm();
  auto cos_phi = (vec1 * vec2) / (d1 * d2);
  if (sanitize_cosine) {
    cos_phi = detail::sanitize_cosine(cos_phi);
  }
  /* force factor */
  auto const fac = forceFactor(cos_phi);
  /* distribute forces */
  auto const v1 = vec1 / d1;
  auto const v2 = vec2 / d2;
  auto f_left = (fac / d1) * (v1 * cos_phi - v2);
  auto f_right = (fac / d2) * (v2 * cos_phi - v1);
  auto f_mid = -(f_left + f_right);
  return std::make_tuple(f_mid, f_left, f_right);
}

#endif /* ANGLE_COMMON_H */
