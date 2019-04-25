/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
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

#include "Rhomboid.hpp"

#include <cmath>

using namespace std;

namespace Shapes {
void Rhomboid::calculate_dist(const Utils::Vector3d &pos, double *dist,
                              double *vec) const {
  double axb[3], bxc[3], axc[3];
  double A, B, C;
  double a_dot_bxc, b_dot_axc, c_dot_axb;
  double tmp;
  double d;

  // calculate a couple of vectors and scalars that are going to be used
  // frequently

  axb[0] = m_a[1] * m_b[2] - m_a[2] * m_b[1];
  axb[1] = m_a[2] * m_b[0] - m_a[0] * m_b[2];
  axb[2] = m_a[0] * m_b[1] - m_a[1] * m_b[0];

  bxc[0] = m_b[1] * m_c[2] - m_b[2] * m_c[1];
  bxc[1] = m_b[2] * m_c[0] - m_b[0] * m_c[2];
  bxc[2] = m_b[0] * m_c[1] - m_b[1] * m_c[0];

  axc[0] = m_a[1] * m_c[2] - m_a[2] * m_c[1];
  axc[1] = m_a[2] * m_c[0] - m_a[0] * m_c[2];
  axc[2] = m_a[0] * m_c[1] - m_a[1] * m_c[0];

  a_dot_bxc = m_a[0] * bxc[0] + m_a[1] * bxc[1] + m_a[2] * bxc[2];
  b_dot_axc = m_b[0] * axc[0] + m_b[1] * axc[1] + m_b[2] * axc[2];
  c_dot_axb = m_c[0] * axb[0] + m_c[1] * axb[1] + m_c[2] * axb[2];

  // represent the distance from pos to ppos as a linear combination of the edge
  // vectors.

  A = (pos[0] - m_pos[0]) * bxc[0] + (pos[1] - m_pos[1]) * bxc[1] +
      (pos[2] - m_pos[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0]) * axc[0] + (pos[1] - m_pos[1]) * axc[1] +
      (pos[2] - m_pos[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0]) * axb[0] + (pos[1] - m_pos[1]) * axb[1] +
      (pos[2] - m_pos[2]) * axb[2];
  C /= c_dot_axb;

  // the coefficients tell whether ppos lies within the cone defined by pos and
  // the adjacent edges

  if (A <= 0 && B <= 0 && C <= 0) {
    vec[0] = pos[0] - m_pos[0];
    vec[1] = pos[1] - m_pos[1];
    vec[2] = pos[2] - m_pos[2];

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for cone at pos+a

  A = (pos[0] - m_pos[0] - m_a[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_a[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_a[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && B <= 0 && C <= 0) {
    vec[0] = pos[0] - m_pos[0] - m_a[0];
    vec[1] = pos[1] - m_pos[1] - m_a[1];
    vec[2] = pos[2] - m_pos[2] - m_a[2];

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for cone at pos+b

  A = (pos[0] - m_pos[0] - m_b[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_b[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_b[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_b[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_b[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_b[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_b[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_b[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_b[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && B >= 0 && C <= 0) {
    vec[0] = pos[0] - m_pos[0] - m_b[0];
    vec[1] = pos[1] - m_pos[1] - m_b[1];
    vec[2] = pos[2] - m_pos[2] - m_b[2];

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for cone at pos+c

  A = (pos[0] - m_pos[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_c[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_c[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && B <= 0 && C >= 0) {
    vec[0] = pos[0] - m_pos[0] - m_c[0];
    vec[1] = pos[1] - m_pos[1] - m_c[1];
    vec[2] = pos[2] - m_pos[2] - m_c[2];

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for cone at m_pos+a+b

  A = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && B >= 0 && C <= 0) {
    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_b[0];
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_b[1];
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_b[2];

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for cone at m_pos+a+c

  A = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && B <= 0 && C >= 0) {
    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_c[0];
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_c[1];
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_c[2];

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for cone at m_pos+a+c

  A = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && B >= 0 && C >= 0) {
    vec[0] = pos[0] - m_pos[0] - m_b[0] - m_c[0];
    vec[1] = pos[1] - m_pos[1] - m_b[1] - m_c[1];
    vec[2] = pos[2] - m_pos[2] - m_b[2] - m_c[2];

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for cone at m_pos+a+b+c

  A = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && B >= 0 && C >= 0) {
    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0];
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1];
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2];

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos, a

  B = (pos[0] - m_pos[0]) * axc[0] + (pos[1] - m_pos[1]) * axc[1] +
      (pos[2] - m_pos[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0]) * axb[0] + (pos[1] - m_pos[1]) * axb[1] +
      (pos[2] - m_pos[2]) * axb[2];
  C /= c_dot_axb;

  if (B <= 0 && C <= 0) {
    tmp = (pos[0] - m_pos[0]) * m_a[0] + (pos[1] - m_pos[1]) * m_a[1] +
          (pos[2] - m_pos[2]) * m_a[2];
    tmp /= m_a[0] * m_a[0] + m_a[1] * m_a[1] + m_a[2] * m_a[2];

    vec[0] = pos[0] - m_pos[0] - m_a[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_a[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_a[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos, b

  A = (pos[0] - m_pos[0]) * bxc[0] + (pos[1] - m_pos[1]) * bxc[1] +
      (pos[2] - m_pos[2]) * bxc[2];
  A /= a_dot_bxc;
  C = (pos[0] - m_pos[0]) * axb[0] + (pos[1] - m_pos[1]) * axb[1] +
      (pos[2] - m_pos[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && C <= 0) {
    tmp = (pos[0] - m_pos[0]) * m_b[0] + (pos[1] - m_pos[1]) * m_b[1] +
          (pos[2] - m_pos[2]) * m_b[2];
    tmp /= m_b[0] * m_b[0] + m_b[1] * m_b[1] + m_b[2] * m_b[2];

    vec[0] = pos[0] - m_pos[0] - m_b[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_b[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_b[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos, c

  A = (pos[0] - m_pos[0]) * bxc[0] + (pos[1] - m_pos[1]) * bxc[1] +
      (pos[2] - m_pos[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0]) * axc[0] + (pos[1] - m_pos[1]) * axc[1] +
      (pos[2] - m_pos[2]) * axc[2];
  B /= b_dot_axc;

  if (A <= 0 && B <= 0) {
    tmp = (pos[0] - m_pos[0]) * m_c[0] + (pos[1] - m_pos[1]) * m_c[1] +
          (pos[2] - m_pos[2]) * m_c[2];
    tmp /= m_c[0] * m_c[0] + m_c[1] * m_c[1] + m_c[2] * m_c[2];

    vec[0] = pos[0] - m_pos[0] - m_c[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_c[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_c[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+a, b

  A = (pos[0] - m_pos[0] - m_a[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2]) * bxc[2];
  A /= a_dot_bxc;
  C = (pos[0] - m_pos[0] - m_a[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && C <= 0) {
    tmp = (pos[0] - m_pos[0] - m_a[0]) * m_b[0] +
          (pos[1] - m_pos[1] - m_a[1]) * m_b[1] +
          (pos[2] - m_pos[2] - m_a[2]) * m_b[2];
    tmp /= m_b[0] * m_b[0] + m_b[1] * m_b[1] + m_b[2] * m_b[2];

    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_b[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_b[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_b[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+a, c

  A = (pos[0] - m_pos[0] - m_a[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_a[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2]) * axc[2];
  B /= b_dot_axc;

  if (A >= 0 && B <= 0) {
    tmp = (pos[0] - m_pos[0] - m_a[0]) * m_c[0] +
          (pos[1] - m_pos[1] - m_a[1]) * m_c[1] +
          (pos[2] - m_pos[2] - m_a[2]) * m_c[2];
    tmp /= m_c[0] * m_c[0] + m_c[1] * m_c[1] + m_c[2] * m_c[2];

    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_c[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_c[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_c[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+b+c, c

  A = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * axc[2];
  B /= b_dot_axc;

  if (A <= 0 && B >= 0) {
    tmp = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * m_c[0] +
          (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * m_c[1] +
          (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * m_c[2];
    tmp /= m_c[0] * m_c[0] + m_c[1] * m_c[1] + m_c[2] * m_c[2];

    vec[0] = pos[0] - m_pos[0] - m_b[0] - m_c[0] - m_c[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_b[1] - m_c[1] - m_c[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_b[2] - m_c[2] - m_c[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+b+c, b

  A = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * bxc[2];
  A /= a_dot_bxc;
  C = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && C >= 0) {
    tmp = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * m_b[0] +
          (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * m_b[1] +
          (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * m_b[2];
    tmp /= m_b[0] * m_b[0] + m_b[1] * m_b[1] + m_b[2] * m_b[2];

    vec[0] = pos[0] - m_pos[0] - m_b[0] - m_c[0] - m_b[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_b[1] - m_c[1] - m_b[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_b[2] - m_c[2] - m_b[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+b+c, a

  B = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * axb[2];
  C /= c_dot_axb;

  if (B >= 0 && C >= 0) {
    tmp = (pos[0] - m_pos[0] - m_b[0] - m_c[0]) * m_a[0] +
          (pos[1] - m_pos[1] - m_b[1] - m_c[1]) * m_a[1] +
          (pos[2] - m_pos[2] - m_b[2] - m_c[2]) * m_a[2];
    tmp /= m_a[0] * m_a[0] + m_a[1] * m_a[1] + m_a[2] * m_a[2];

    vec[0] = pos[0] - m_pos[0] - m_b[0] - m_c[0] - m_a[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_b[1] - m_c[1] - m_a[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_b[2] - m_c[2] - m_a[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+a+b, a

  B = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * axb[2];
  C /= c_dot_axb;

  if (B >= 0 && C <= 0) {
    tmp = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * m_a[0] +
          (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * m_a[1] +
          (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * m_a[2];
    tmp /= m_a[0] * m_a[0] + m_a[1] * m_a[1] + m_a[2] * m_a[2];

    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_a[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_a[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_a[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+a+b, c

  A = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * axc[2];
  B /= b_dot_axc;

  if (A >= 0 && B >= 0) {
    tmp = (pos[0] - m_pos[0] - m_a[0] - m_b[0]) * m_c[0] +
          (pos[1] - m_pos[1] - m_a[1] - m_b[1]) * m_c[1] +
          (pos[2] - m_pos[2] - m_a[2] - m_b[2]) * m_c[2];
    tmp /= m_c[0] * m_c[0] + m_c[1] * m_c[1] + m_c[2] * m_c[2];

    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+a+c, a

  B = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * axc[2];
  B /= b_dot_axc;
  C = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * axb[2];
  C /= c_dot_axb;

  if (B <= 0 && C >= 0) {
    tmp = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * m_a[0] +
          (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * m_a[1] +
          (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * m_a[2];
    tmp /= m_a[0] * m_a[0] + m_a[1] * m_a[1] + m_a[2] * m_a[2];

    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_c[0] - m_a[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_c[1] - m_a[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_c[2] - m_a[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for prism at edge m_pos+a+c, b

  A = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * bxc[2];
  A /= a_dot_bxc;
  C = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && C >= 0) {
    tmp = (pos[0] - m_pos[0] - m_a[0] - m_c[0]) * m_b[0] +
          (pos[1] - m_pos[1] - m_a[1] - m_c[1]) * m_b[1] +
          (pos[2] - m_pos[2] - m_a[2] - m_c[2]) * m_b[2];
    tmp /= m_b[0] * m_b[0] + m_b[1] * m_b[1] + m_b[2] * m_b[2];

    vec[0] = pos[0] - m_pos[0] - m_a[0] - m_c[0] - m_b[0] * tmp;
    vec[1] = pos[1] - m_pos[1] - m_a[1] - m_c[1] - m_b[1] * tmp;
    vec[2] = pos[2] - m_pos[2] - m_a[2] - m_c[2] - m_b[2] * tmp;

    *dist =
        m_direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return;
  }

  // check for face with normal -axb

  *dist = (pos[0] - m_pos[0]) * axb[0] + (pos[1] - m_pos[1]) * axb[1] +
          (pos[2] - m_pos[2]) * axb[2];
  if (c_dot_axb > 0.0)
    *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
    *dist /= tmp;

    vec[0] = -*dist * axb[0] / tmp;
    vec[1] = -*dist * axb[1] / tmp;
    vec[2] = -*dist * axb[2] / tmp;
    if (c_dot_axb < 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist *= m_direction;

    return;
  }

  // calculate distance to face with normal axc

  *dist = (pos[0] - m_pos[0]) * axc[0] + (pos[1] - m_pos[1]) * axc[1] +
          (pos[2] - m_pos[2]) * axc[2];
  if (b_dot_axc > 0.0)
    *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(axc[0] * axc[0] + axc[1] * axc[1] + axc[2] * axc[2]);
    *dist /= tmp;

    vec[0] = *dist * axc[0] / tmp;
    vec[1] = *dist * axc[1] / tmp;
    vec[2] = *dist * axc[2] / tmp;
    if (b_dot_axc > 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist *= m_direction;

    return;
  }

  // calculate distance to face with normal -bxc

  *dist = (pos[0] - m_pos[0]) * bxc[0] + (pos[1] - m_pos[1]) * bxc[1] +
          (pos[2] - m_pos[2]) * bxc[2];
  if (a_dot_bxc > 0.0)
    *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(bxc[0] * bxc[0] + bxc[1] * bxc[1] + bxc[2] * bxc[2]);
    *dist /= tmp;

    vec[0] = -*dist * bxc[0] / tmp;
    vec[1] = -*dist * bxc[1] / tmp;
    vec[2] = -*dist * bxc[2] / tmp;
    if (a_dot_bxc < 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist *= m_direction;

    return;
  }

  // calculate distance to face with normal axb

  *dist = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axb[0] +
          (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axb[1] +
          (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axb[2];
  if (c_dot_axb < 0.0)
    *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
    *dist /= tmp;

    vec[0] = *dist * axb[0] / tmp;
    vec[1] = *dist * axb[1] / tmp;
    vec[2] = *dist * axb[2] / tmp;
    if (c_dot_axb < 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist *= m_direction;

    return;
  }

  // calculate distance to face with normal -axc

  *dist = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axc[0] +
          (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axc[1] +
          (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axc[2];
  if (b_dot_axc < 0.0)
    *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(axc[0] * axc[0] + axc[1] * axc[1] + axc[2] * axc[2]);
    *dist /= tmp;

    vec[0] = -*dist * axc[0] / tmp;
    vec[1] = -*dist * axc[1] / tmp;
    vec[2] = -*dist * axc[2] / tmp;
    if (b_dot_axc > 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist *= m_direction;

    return;
  }

  // calculate distance to face with normal bxc

  *dist = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * bxc[0] +
          (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * bxc[1] +
          (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * bxc[2];
  if (a_dot_bxc < 0.0)
    *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(bxc[0] * bxc[0] + bxc[1] * bxc[1] + bxc[2] * bxc[2]);
    *dist /= tmp;

    vec[0] = *dist * bxc[0] / tmp;
    vec[1] = *dist * bxc[1] / tmp;
    vec[2] = *dist * bxc[2] / tmp;
    if (a_dot_bxc < 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist *= m_direction;

    return;
  }

  // ppos lies within rhomboid. Find nearest wall for interaction.

  // check for face with normal -axb

  *dist = (pos[0] - m_pos[0]) * axb[0] + (pos[1] - m_pos[1]) * axb[1] +
          (pos[2] - m_pos[2]) * axb[2];
  if (c_dot_axb > 0.0)
    *dist *= -1.;
  tmp = sqrt(axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
  *dist /= tmp;

  vec[0] = -*dist * axb[0] / tmp;
  vec[1] = -*dist * axb[1] / tmp;
  vec[2] = -*dist * axb[2] / tmp;

  if (c_dot_axb < 0.0) {
    vec[0] *= -1.;
    vec[1] *= -1.;
    vec[2] *= -1.;
  }

  *dist *= m_direction;

  // calculate distance to face with normal axc

  d = (pos[0] - m_pos[0]) * axc[0] + (pos[1] - m_pos[1]) * axc[1] +
      (pos[2] - m_pos[2]) * axc[2];
  if (b_dot_axc > 0.0)
    d *= -1.;
  tmp = sqrt(axc[0] * axc[0] + axc[1] * axc[1] + axc[2] * axc[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = d * axc[0] / tmp;
    vec[1] = d * axc[1] / tmp;
    vec[2] = d * axc[2] / tmp;

    if (b_dot_axc > 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist = m_direction * d;
  }

  // calculate distance to face with normal -bxc

  d = (pos[0] - m_pos[0]) * bxc[0] + (pos[1] - m_pos[1]) * bxc[1] +
      (pos[2] - m_pos[2]) * bxc[2];
  if (a_dot_bxc > 0.0)
    d *= -1.;
  tmp = sqrt(bxc[0] * bxc[0] + bxc[1] * bxc[1] + bxc[2] * bxc[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = -d * bxc[0] / tmp;
    vec[1] = -d * bxc[1] / tmp;
    vec[2] = -d * bxc[2] / tmp;

    if (a_dot_bxc < 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist = m_direction * d;
  }

  // calculate distance to face with normal axb

  d = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axb[2];
  if (c_dot_axb < 0.0)
    d *= -1.;
  tmp = sqrt(axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = d * axb[0] / tmp;
    vec[1] = d * axb[1] / tmp;
    vec[2] = d * axb[2] / tmp;

    if (c_dot_axb < 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist = m_direction * d;
  }

  // calculate distance to face with normal -axc

  d = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axc[2];
  if (b_dot_axc < 0.0)
    d *= -1.;
  tmp = sqrt(axc[0] * axc[0] + axc[1] * axc[1] + axc[2] * axc[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = -d * axc[0] / tmp;
    vec[1] = -d * axc[1] / tmp;
    vec[2] = -d * axc[2] / tmp;

    if (b_dot_axc > 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist = m_direction * d;
  }

  // calculate distance to face with normal bxc

  d = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * bxc[2];
  if (a_dot_bxc < 0.0)
    d *= -1.;
  tmp = sqrt(bxc[0] * bxc[0] + bxc[1] * bxc[1] + bxc[2] * bxc[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = d * bxc[0] / tmp;
    vec[1] = d * bxc[1] / tmp;
    vec[2] = d * bxc[2] / tmp;

    if (a_dot_bxc < 0.0) {
      vec[0] *= -1.;
      vec[1] *= -1.;
      vec[2] *= -1.;
    }

    *dist = m_direction * d;
  }
}

} // namespace Shapes
