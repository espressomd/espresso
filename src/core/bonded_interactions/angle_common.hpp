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
#include "particle_data.hpp"
#include "utils.hpp"

/** Compute a three-body angle interaction force.
 *  @param[in]  p_mid        Second/middle particle.
 *  @param[in]  p_left       First/left particle.
 *  @param[in]  p_right      Third/right particle.
 *  @param[in]  forceFactor  Angle bending constant.
 *  @param[out] force1       Force on particle 1.
 *  @param[out] force2       Force on particle 2.
 */
template <typename ForceFactor>
void calc_angle_generic_force(Particle const *p_mid, Particle const *p_left,
                              Particle const *p_right, ForceFactor forceFactor,
                              double force1[3], double force2[3]) {
  /* vector from p_left to p_mid */
  auto vec1 = get_mi_vector(p_mid->r.p, p_left->r.p);
  auto d1i = 1.0 / vec1.norm();
  vec1 *= d1i;
  /* vector from p_mid to p_right */
  auto vec2 = get_mi_vector(p_right->r.p, p_mid->r.p);
  auto d2i = 1.0 / vec2.norm();
  vec2 *= d2i;
  /* scalar product of vec1 and vec2 */
  auto cosine = scalar(vec1, vec2);
  /* force factor */
  auto fac = forceFactor(cosine);
  /* force calculation */
  auto f1 = fac * (cosine * vec1 - vec2) * d1i;
  auto f2 = fac * (cosine * vec2 - vec1) * d2i;
  for (int j = 0; j < 3; j++) {
    force1[j] = (f1 - f2)[j];
    force2[j] = -f1[j];
  }
}

/** Compute the forces of a three-body bonded potential.
 *  @param[in]  p_mid        Second/middle particle.
 *  @param[in]  p_left       First/left particle.
 *  @param[in]  p_right      Third/right particle.
 *  @param[in]  forceFactor  Angle bending constant.
 *  @param[out] force1       Force on particle 1.
 *  @param[out] force2       Force on particle 2.
 *  @param[out] force3       Force on particle 3.
 */
template <typename ForceFactor>
void calc_angle_generic_3body_forces(Particle const *p_mid,
                                     Particle const *p_left,
                                     Particle const *p_right,
                                     ForceFactor forceFactor, Vector3d &force1,
                                     Vector3d &force2, Vector3d &force3) {
  auto vec21 = -get_mi_vector(p_mid->r.p, p_left->r.p);
  auto vec31 = get_mi_vector(p_right->r.p, p_mid->r.p);
  auto vec21_sqr = vec21.norm2();
  auto vec31_sqr = vec31.norm2();
  auto vec21_magn = sqrt(vec21_sqr);
  auto vec31_magn = sqrt(vec31_sqr);
  auto cos_phi = scalar(vec21, vec31) / (vec21_magn * vec31_magn);
  auto sin_phi = sqrt(1.0 - Utils::sqr(cos_phi));
  /* force factor */
  auto fac = forceFactor(cos_phi, sin_phi);
  /* force calculation */
  auto fj = vec31 / (vec21_magn * vec31_magn) - cos_phi * vec21 / vec21_sqr;
  auto fk = vec21 / (vec21_magn * vec31_magn) - cos_phi * vec31 / vec31_sqr;
  // note that F1 = -(F2 + F3) in analytical case
  force1 = -fac * (fj + fk);
  force2 = fac * fj;
  force3 = fac * fk;
}

#endif /* ANGLE_COMMON_H */
