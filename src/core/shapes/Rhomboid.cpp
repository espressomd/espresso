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
void Rhomboid::calculate_dist(const Utils::Vector3d &pos, double &dist,
                              Utils::Vector3d &vec) const {
  double tmp;
  double d;

  // calculate a couple of vectors and scalars that are going to be used
  // frequently

  auto const axb = vector_product(m_a, m_b);
  auto const bxc = vector_product(m_b, m_c);
  auto const axc = vector_product(m_a, m_c);

  auto const a_dot_bxc = m_a * bxc;
  auto const b_dot_axc = m_b * axc;
  auto const c_dot_axb = m_c * axb;

  // represent the distance from pos to ppos as a linear combination of the edge
  // vectors.

  auto const dpos = pos - m_pos;
  auto A = dpos * bxc / a_dot_bxc;
  auto B = dpos * axc / b_dot_axc;
  auto C = dpos * axb / c_dot_axb;

  // the coefficients tell whether ppos lies within the cone defined by pos and
  // the adjacent edges

  if (A <= 0 && B <= 0 && C <= 0) {
    vec = dpos;
    dist = m_direction * vec.norm();
    return;
  }

  // check for cone at pos+a

  auto const dpos_a = dpos - m_a;
  A = dpos_a * bxc / a_dot_bxc;
  B = dpos_a * axc / b_dot_axc;
  C = dpos_a * axb / c_dot_axb;

  if (A >= 0 && B <= 0 && C <= 0) {
    vec = dpos_a;
    dist = m_direction * vec.norm();
    return;
  }

  // check for cone at pos+b

  auto const dpos_b = dpos - m_b;
  A = dpos_b * bxc / a_dot_bxc;
  B = dpos_b * axc / b_dot_axc;
  C = dpos_b * axb / c_dot_axb;

  if (A <= 0 && B >= 0 && C <= 0) {
    vec = dpos_b;
    dist = m_direction * vec.norm();
    return;
  }

  // check for cone at pos+c

  auto const dpos_c = dpos - m_c;
  A = dpos_c * bxc / a_dot_bxc;
  B = dpos_c * axc / b_dot_axc;
  C = dpos_c * axb / c_dot_axb;

  if (A <= 0 && B <= 0 && C >= 0) {
    vec = dpos_c;
    dist = m_direction * vec.norm();
    return;
  }

  // check for cone at m_pos+a+b

  auto const dpos_ab = dpos - m_a - m_b;
  A = dpos_ab * bxc / a_dot_bxc;
  B = dpos_ab * axc / b_dot_axc;
  C = dpos_ab * axb / c_dot_axb;

  if (A >= 0 && B >= 0 && C <= 0) {
    vec = dpos_ab;
    dist = m_direction * vec.norm();
    return;
  }

  // check for cone at m_pos+a+c

  auto const dpos_ac = dpos - m_a - m_c;
  A = dpos_ac * bxc / a_dot_bxc;
  B = dpos_ac * axc / b_dot_axc;
  C = dpos_ac * axb / c_dot_axb;

  if (A >= 0 && B <= 0 && C >= 0) {
    vec = dpos_ac;
    dist = m_direction * vec.norm();
    return;
  }

  // check for cone at m_pos+b+c

  auto const dpos_bc = dpos - m_b - m_c;
  A = dpos_bc * bxc / a_dot_bxc;
  B = dpos_bc * axc / b_dot_axc;
  C = dpos_bc * axb / c_dot_axb;

  if (A <= 0 && B >= 0 && C >= 0) {
    vec = dpos_bc;
    dist = m_direction * vec.norm();
    return;
  }

  // check for cone at m_pos+a+b+c

  auto const dpos_abc = dpos - m_a - m_b - m_c;
  A = dpos_abc * bxc / a_dot_bxc;
  B = dpos_abc * axc / b_dot_axc;
  C = dpos_abc * axb / c_dot_axb;

  if (A >= 0 && B >= 0 && C >= 0) {
    vec = dpos_abc;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_a * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_b * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_c * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_a - m_b * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_a - m_c * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_b - m_c - m_c * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_b - m_c - m_b * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_b - m_c - m_a * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_a - m_b - m_a * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_a - m_b - m_c * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_a - m_c - m_a * tmp;
    dist = m_direction * vec.norm();
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

    vec = pos - m_pos - m_a - m_c - m_b * tmp;
    dist = m_direction * vec.norm();
    return;
  }

  // check for face with normal -axb

  dist = (pos[0] - m_pos[0]) * axb[0] + (pos[1] - m_pos[1]) * axb[1] +
         (pos[2] - m_pos[2]) * axb[2];
  if (c_dot_axb > 0.0) {
    dist *= -1.;
  }

  if (dist >= 0) {
    tmp = axb.norm();
    dist /= tmp;

    vec = (-dist / tmp) * axb;
    if (c_dot_axb < 0.0) {
      vec *= -1.;
    }

    dist *= m_direction;

    return;
  }

  // calculate distance to face with normal axc

  dist = (pos[0] - m_pos[0]) * axc[0] + (pos[1] - m_pos[1]) * axc[1] +
         (pos[2] - m_pos[2]) * axc[2];
  if (b_dot_axc > 0.0) {
    dist *= -1.;
  }

  if (dist >= 0) {
    tmp = axc.norm();
    dist /= tmp;

    vec = (dist / tmp) * axc;
    if (b_dot_axc > 0.0) {
      vec *= -1.;
    }

    dist *= m_direction;

    return;
  }

  // calculate distance to face with normal -bxc

  dist = (pos[0] - m_pos[0]) * bxc[0] + (pos[1] - m_pos[1]) * bxc[1] +
         (pos[2] - m_pos[2]) * bxc[2];
  if (a_dot_bxc > 0.0) {
    dist *= -1.;
  }

  if (dist >= 0) {
    tmp = bxc.norm();
    dist /= tmp;

    vec = -(dist / tmp) * bxc;
    if (a_dot_bxc < 0.0) {
      vec *= -1.;
    }

    dist *= m_direction;

    return;
  }

  // calculate distance to face with normal axb

  dist = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axb[0] +
         (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axb[1] +
         (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axb[2];
  if (c_dot_axb < 0.0) {
    dist *= -1.;
  }

  if (dist >= 0) {
    tmp = axb.norm();
    dist /= tmp;

    vec = (dist / tmp) * axb;
    if (c_dot_axb < 0.0) {
      vec *= -1.;
    }

    dist *= m_direction;

    return;
  }

  // calculate distance to face with normal -axc

  dist = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axc[0] +
         (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axc[1] +
         (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axc[2];
  if (b_dot_axc < 0.0) {
    dist *= -1.;
  }

  if (dist >= 0) {
    tmp = axc.norm();
    dist /= tmp;

    vec = -(dist / tmp) * axc;
    if (b_dot_axc > 0.0) {
      vec *= -1.;
    }

    dist *= m_direction;

    return;
  }

  // calculate distance to face with normal bxc

  dist = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * bxc[0] +
         (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * bxc[1] +
         (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * bxc[2];
  if (a_dot_bxc < 0.0) {
    dist *= -1.;
  }

  if (dist >= 0) {
    tmp = bxc.norm();
    dist /= tmp;

    vec = (dist / tmp) * bxc;
    if (a_dot_bxc < 0.0) {
      vec *= -1.;
    }

    dist *= m_direction;

    return;
  }

  // ppos lies within rhomboid. Find nearest wall for interaction.

  // check for face with normal -axb

  dist = (pos[0] - m_pos[0]) * axb[0] + (pos[1] - m_pos[1]) * axb[1] +
         (pos[2] - m_pos[2]) * axb[2];
  if (c_dot_axb > 0.0) {
    dist *= -1.;
  }
  tmp = axb.norm();
  dist /= tmp;

  vec = -(dist / tmp) * axb;

  if (c_dot_axb < 0.0) {
    vec *= -1.;
  }

  dist *= m_direction;

  // calculate distance to face with normal axc

  d = (pos[0] - m_pos[0]) * axc[0] + (pos[1] - m_pos[1]) * axc[1] +
      (pos[2] - m_pos[2]) * axc[2];
  if (b_dot_axc > 0.0) {
    d *= -1.;
  }
  tmp = axc.norm();
  d /= tmp;

  if (abs(d) < abs(dist)) {
    vec = (d / tmp) * axc;

    if (b_dot_axc > 0.0) {
      vec *= -1.;
    }

    dist = m_direction * d;
  }

  // calculate distance to face with normal -bxc

  d = (pos[0] - m_pos[0]) * bxc[0] + (pos[1] - m_pos[1]) * bxc[1] +
      (pos[2] - m_pos[2]) * bxc[2];
  if (a_dot_bxc > 0.0)
    d *= -1.;
  tmp = bxc.norm();
  d /= tmp;

  if (abs(d) < abs(dist)) {
    vec = -(d / tmp) * bxc;

    if (a_dot_bxc < 0.0) {
      vec *= -1.;
    }

    dist = m_direction * d;
  }

  // calculate distance to face with normal axb

  d = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axb[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axb[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axb[2];
  if (c_dot_axb < 0.0) {
    d *= -1.;
  }
  tmp = axb.norm();
  d /= tmp;

  if (abs(d) < abs(dist)) {
    vec = (d / tmp) * axb;

    if (c_dot_axb < 0.0) {
      vec *= -1.;
    }

    dist = m_direction * d;
  }

  // calculate distance to face with normal -axc

  d = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * axc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * axc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * axc[2];
  if (b_dot_axc < 0.0)
    d *= -1.;
  tmp = axc.norm();
  d /= tmp;

  if (abs(d) < abs(dist)) {
    vec = -(d / tmp) * axc;

    if (b_dot_axc > 0.0) {
      vec *= -1.;
    }

    dist = m_direction * d;
  }

  // calculate distance to face with normal bxc

  d = (pos[0] - m_pos[0] - m_a[0] - m_b[0] - m_c[0]) * bxc[0] +
      (pos[1] - m_pos[1] - m_a[1] - m_b[1] - m_c[1]) * bxc[1] +
      (pos[2] - m_pos[2] - m_a[2] - m_b[2] - m_c[2]) * bxc[2];
  if (a_dot_bxc < 0.0)
    d *= -1.;
  tmp = bxc.norm();
  d /= tmp;

  if (abs(d) < abs(dist)) {
    vec = (d / tmp) * bxc;

    if (a_dot_bxc < 0.0) {
      vec *= -1.;
    }

    dist = m_direction * d;
  }
}

} // namespace Shapes
