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
 *  @f$ \vec{r_{ji}} @f$ (from particle @f$ i @f$ to particle @f$ j @f$)
 *  and @f$ \vec{r_{kj}} @f$, and their normalization constants.
 *
 *  @param[in]  r_mid            Position of second/middle particle.
 *  @param[in]  r_left           Position of first/left particle.
 *  @param[in]  r_right          Position of third/right particle.
 *  @param[in]  sanitize_cosine  Sanitize the cosine of the angle.
 *  @return @f$ \vec{r_{ji}} @f$, @f$ \vec{r_{kj}} @f$,
 *          @f$ \left\|\vec{r_{ji}}\right\|^{-1} @f$,
 *          @f$ \left\|\vec{r_{kj}}\right\|^{-1} @f$,
 *          @f$ \cos(\theta_{ijk}) @f$
 */
inline std::tuple<Vector3d, Vector3d, double, double, double>
calc_vectors_and_cosine(Vector3d const &r_mid, Vector3d const &r_left,
                        Vector3d const &r_right, bool sanitize_cosine = false) {
  /* vector from p_left to p_mid */
  auto vec1 = get_mi_vector(r_mid, r_left);
  auto const d1i = 1.0 / vec1.norm();
  vec1 *= d1i;
  /* vector from p_mid to p_right */
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
 *  @param[out] f_mid            Force on the second/middle particle.
 *  @param[out] f_left           Force on the first/left particle.
 *  tparam      ForceFactor      Function evaluating the angle force term
 *                               for a given angle.
 */
template <typename ForceFactor>
void calc_angle_generic_force(Vector3d const &r_mid, Vector3d const &r_left,
                              Vector3d const &r_right, ForceFactor forceFactor,
                              double f_mid[3], double f_left[3],
                              bool sanitize_cosine) {
  Vector3d vec1, vec2;
  double d1i, d2i, cosine;
  std::tie(vec1, vec2, d1i, d2i, cosine) =
      calc_vectors_and_cosine(r_mid, r_left, r_right, sanitize_cosine);
  /* force factor */
  auto const fac = forceFactor(cosine);
  /* force calculation on particles 1 and 3 */
  auto const f1 = fac * (cosine * vec1 - vec2) * d1i;
  auto const f3 = fac * (cosine * vec2 - vec1) * d2i;
  for (int j = 0; j < 3; j++) {
    f_mid[j] = (f1 - f3)[j];
    f_left[j] = -f1[j];
  }
}

/** Compute the forces of a three-body bonded potential.
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
std::tuple<Vector3d, Vector3d, Vector3d>
calc_angle_generic_3body_forces(Vector3d const &r_mid, Vector3d const &r_left,
                                Vector3d const &r_right,
                                ForceFactor forceFactor, bool sanitize_cosine) {
  auto const vec21 = get_mi_vector(r_left, r_mid);
  auto const vec31 = get_mi_vector(r_right, r_mid);
  auto const vec21_sqr = vec21.norm2();
  auto const vec31_sqr = vec31.norm2();
  auto const vec21_len = sqrt(vec21_sqr);
  auto const vec31_len = sqrt(vec31_sqr);
  auto cos_phi = (vec21 * vec31) / (vec21_len * vec31_len);
  if (sanitize_cosine) {
    if (cos_phi > TINY_COS_VALUE)
      cos_phi = TINY_COS_VALUE;
    if (cos_phi < -TINY_COS_VALUE)
      cos_phi = -TINY_COS_VALUE;
  }
  /* force factor */
  auto const fac = forceFactor(cos_phi);
  /* force calculation */
  auto const fj = vec31 / (vec21_len * vec31_len) - cos_phi * vec21 / vec21_sqr;
  auto const fk = vec21 / (vec21_len * vec31_len) - cos_phi * vec31 / vec31_sqr;
  // note that F1 = -(F2 + F3) in analytical case
  auto force_mid = -fac * (fj + fk);
  auto force_left = fac * fj;
  auto force_right = fac * fk;
  return std::make_tuple(force_mid, force_left, force_right);
}

#endif /* ANGLE_COMMON_H */
