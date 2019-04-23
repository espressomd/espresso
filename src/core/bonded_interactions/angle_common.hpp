/*
  Copyright (C) 2010-2019 The ESPResSo project
  Copyright (C) 2002-2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ANGLE_COMMON_H
#define ANGLE_COMMON_H
/** \file
 *  Common code for functions calculating angle forces.
 */

#include "grid.hpp"
#include <tuple>

/** Compute the cosine of the angle between three particles.
 *
 *  Also return all intermediate quantities: normalized vectors
 *  @f$ \vec{r_{ij}} @f$ (from particle @f$ j @f$ to particle @f$ i @f$)
 *  and @f$ \vec{r_{kj}} @f$, and their normalization constants.
 *
 *  @param[in]  r_mid            Position of second/middle particle.
 *  @param[in]  r_left           Position of first/left particle.
 *  @param[in]  r_right          Position of third/right particle.
 *  @param[in]  sanitize_cosine  Sanitize the cosine of the angle.
 *  @return @f$ \vec{r_{ij}} @f$, @f$ \vec{r_{kj}} @f$,
 *          @f$ \left\|\vec{r_{ij}}\right\|^{-1} @f$,
 *          @f$ \left\|\vec{r_{kj}}\right\|^{-1} @f$,
 *          @f$ \cos(\theta_{ijk}) @f$
 */
inline std::tuple<Utils::Vector3d, Utils::Vector3d, double, double, double>
calc_vectors_and_cosine(Utils::Vector3d const &r_mid,
                        Utils::Vector3d const &r_left,
                        Utils::Vector3d const &r_right,
                        bool sanitize_cosine = false) {
  /* normalized vector from p_mid to p_left */
  auto vec1 = get_mi_vector(r_left, r_mid);
  auto const d1i = 1.0 / vec1.norm();
  vec1 *= d1i;
  /* normalized vector from p_mid to p_right */
  auto vec2 = get_mi_vector(r_right, r_mid);
  auto const d2i = 1.0 / vec2.norm();
  vec2 *= d2i;
  /* cosine of the angle between vec1 and vec2 */
  auto cosine = vec1 * vec2;
  if (sanitize_cosine) {
    if (cosine > TINY_COS_VALUE)
      cosine = TINY_COS_VALUE;
    if (cosine < -TINY_COS_VALUE)
      cosine = -TINY_COS_VALUE;
  }
  return std::make_tuple(vec1, vec2, d1i, d2i, cosine);
}

/** Compute a three-body angle interaction force.
 *
 *  See the details in @ref bondedIA_angle_force. The @f$ K(\theta_{ijk}) @f$
 *  term is provided as a lambda function in @p forceFactor.
 *
 *  @param[in]  r_mid            Position of second/middle particle.
 *  @param[in]  r_left           Position of first/left particle.
 *  @param[in]  r_right          Position of third/right particle.
 *  @param[in]  forceFactor      Angle force term.
 *  @param[in]  sanitize_cosine  Sanitize the cosine of the angle.
 *  tparam      ForceFactor      Function evaluating the angle force term
 *                               for a given angle.
 *  @return Forces on the second, first and third particles, in that order.
 */
template <typename ForceFactor>
std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
calc_angle_generic_force(Utils::Vector3d const &r_mid,
                         Utils::Vector3d const &r_left,
                         Utils::Vector3d const &r_right,
                         ForceFactor forceFactor, bool sanitize_cosine) {
  Utils::Vector3d vec1, vec2;
  double d1i, d2i, cosine;
  std::tie(vec1, vec2, d1i, d2i, cosine) =
      calc_vectors_and_cosine(r_mid, r_left, r_right, sanitize_cosine);
  /* force factor */
  auto const fac = forceFactor(cosine);
  /* distribute forces */
  auto f_left = (fac * d1i) * (vec1 * cosine - vec2);
  auto f_right = (fac * d2i) * (vec2 * cosine - vec1);
  auto f_mid = -(f_left + f_right);
  return std::make_tuple(f_mid, f_left, f_right);
}

#endif /* ANGLE_COMMON_H */
