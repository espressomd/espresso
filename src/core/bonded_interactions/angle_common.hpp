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

#include "utils.hpp"
#include "grid.hpp"

inline void calc_angle_vector(const Vector3d &r_mid, const Vector3d &r_out,
                              double vec[3], double &di) {
  /* vector from r_out to r_mid */
  get_mi_vector(vec, r_mid, r_out);
  di = 1.0 / sqrt(sqrlen(vec));
  for (int j = 0; j < 3; j++)
    vec[j] *= di;
}

inline void calc_angle_force(double force1[3], double force2[3],
                             const double vec1[3], const double vec2[3],
                             double d1i, double d2i,
                             double cosine, double fac) {
  for (int j = 0; j < 3; j++) {
    double f1 = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    double f2 = fac * (cosine * vec2[j] - vec1[j]) * d2i;
    force1[j] = f1 - f2;
    force2[j] = -f1;
  }
}

inline void calc_angle_3body_vector(const Vector3d &r_mid,
                                    const Vector3d &r_left,
                                    const Vector3d &r_right,
                                    double &cos_phi,
                                    double &sin_phi,
                                    double vec21[3],
                                    double vec31[3],
                                    double &vec21_sqr,
                                    double &vec31_sqr,
                                    double &vec21_magn,
                                    double &vec31_magn) {
  get_mi_vector(vec21, r_mid, r_left);
  for (int j = 0; j < 3; j++)
    vec21[j] = -vec21[j];
  get_mi_vector(vec31, r_right, r_mid);
  vec21_sqr = sqrlen(vec21);
  vec21_magn = sqrt(vec21_sqr);
  vec31_sqr = sqrlen(vec31);
  vec31_magn = sqrt(vec31_sqr);
  cos_phi = scalar(vec21, vec31) / (vec21_magn * vec31_magn);
  sin_phi = sqrt(1.0 - Utils::sqr(cos_phi));
}

inline void calc_angle_3body_force(double cos_phi,
                                   double fac,
                                   const double vec21[3],
                                   const double vec31[3],
                                   double vec21_sqr,
                                   double vec31_sqr,
                                   double vec21_magn,
                                   double vec31_magn,
                                   double force1[3],
                                   double force2[3],
                                   double force3[3]) {
  double fj, fk;
  for (int j = 0; j < 3; j++) {
    fj = vec31[j] / (vec21_magn * vec31_magn) - cos_phi * vec21[j] / vec21_sqr;
    fk = vec21[j] / (vec21_magn * vec31_magn) - cos_phi * vec31[j] / vec31_sqr;

    // note that F1 = -(F2 + F3) in analytical case
    force1[j] -= fac * (fj + fk);
    force2[j] += fac * fj;
    force3[j] += fac * fk;
  }
}

#endif /* ANGLE_COMMON_H */
