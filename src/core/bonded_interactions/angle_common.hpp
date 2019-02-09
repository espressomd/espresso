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
 *  @param[in]  r_mid            Position of second/middle particle.
 *  @param[in]  r_left           Position of first/left particle.
 *  @param[in]  r_right          Position of third/right particle.
 *  @param[in]  forceFactor      Angle bending constant.
 *  @param[in]  sanitize_cosine  Sanitize the cosine of the angle.
 *  @param[out] force1           Force on particle 1.
 *  @param[out] force2           Force on particle 2.
 */
template <typename ForceFactor>
void calc_angle_generic_force(Vector3d const &r_mid, Vector3d const &r_left,
                              Vector3d const &r_right, ForceFactor forceFactor,
                              double force1[3], double force2[3],
                              bool sanitize_cosine = false) {
  Vector3d vec1, vec2;
  double d1i, d2i, cosine;
  std::tie(vec1, vec2, d1i, d2i, cosine) =
      calc_vectors_and_cosine(r_mid, r_left, r_right, sanitize_cosine);
  /* force factor */
  auto const fac = forceFactor(cosine);
  /* force calculation */
  auto const f1 = fac * (cosine * vec1 - vec2) * d1i;
  auto const f2 = fac * (cosine * vec2 - vec1) * d2i;
  for (int j = 0; j < 3; j++) {
    force1[j] = (f1 - f2)[j];
    force2[j] = -f1[j];
  }
}

/** Compute the forces of a three-body bonded potential.
 *  @param[in]  r_mid        Position of second/middle particle.
 *  @param[in]  r_left       Position of first/left particle.
 *  @param[in]  r_right      Position of third/right particle.
 *  @param[in]  forceFactor  Angle bending constant.
 *  @return Forces on particles 1, 2 and 3.
 */
template <typename ForceFactor>
std::tuple<Vector3d, Vector3d, Vector3d>
calc_angle_generic_3body_forces(Vector3d const &r_mid, Vector3d const &r_left,
                                Vector3d const &r_right,
                                ForceFactor forceFactor) {
  auto const vec21 = get_mi_vector(r_left, r_mid);
  auto const vec31 = get_mi_vector(r_right, r_mid);
  auto const vec21_sqr = vec21.norm2();
  auto const vec31_sqr = vec31.norm2();
  auto const vec21_len = sqrt(vec21_sqr);
  auto const vec31_len = sqrt(vec31_sqr);
  auto cos_phi = (vec21 * vec31) / (vec21_len * vec31_len);
  auto sin_phi = sqrt(1.0 - Utils::sqr(cos_phi));
  /* force factor */
  auto const fac = forceFactor(cos_phi, sin_phi);
  /* force calculation */
  auto const fj = vec31 / (vec21_len * vec31_len) - cos_phi * vec21 / vec21_sqr;
  auto const fk = vec21 / (vec21_len * vec31_len) - cos_phi * vec31 / vec31_sqr;
  // note that F1 = -(F2 + F3) in analytical case
  auto force1 = -fac * (fj + fk);
  auto force2 = fac * fj;
  auto force3 = fac * fk;
  return std::make_tuple(force1, force2, force3);
}

#endif /* ANGLE_COMMON_H */
