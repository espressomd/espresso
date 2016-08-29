/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
int Rhomboid::calculate_dist(const double *ppos, double *dist, double *vec) const {
  double axb[3], bxc[3], axc[3];
  double A, B, C;
  double a_dot_bxc, b_dot_axc, c_dot_axb;
  double tmp;
  double d;

  // calculate a couple of vectors and scalars that are going to be used
  // frequently

  axb[0] = a[1] * b[2] - a[2] * b[1];
  axb[1] = a[2] * b[0] - a[0] * b[2];
  axb[2] = a[0] * b[1] - a[1] * b[0];

  bxc[0] = b[1] * c[2] - b[2] * c[1];
  bxc[1] = b[2] * c[0] - b[0] * c[2];
  bxc[2] = b[0] * c[1] - b[1] * c[0];

  axc[0] = a[1] * c[2] - a[2] * c[1];
  axc[1] = a[2] * c[0] - a[0] * c[2];
  axc[2] = a[0] * c[1] - a[1] * c[0];

  a_dot_bxc = a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2];
  b_dot_axc = b[0] * axc[0] + b[1] * axc[1] + b[2] * axc[2];
  c_dot_axb = c[0] * axb[0] + c[1] * axb[1] + c[2] * axb[2];

  // represent the distance from pos to ppos as a linear combination of the edge
  // vectors.

  A = (ppos[0] - pos[0]) * bxc[0] + (ppos[1] - pos[1]) * bxc[1] +
      (ppos[2] - pos[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0]) * axc[0] + (ppos[1] - pos[1]) * axc[1] +
      (ppos[2] - pos[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0]) * axb[0] + (ppos[1] - pos[1]) * axb[1] +
      (ppos[2] - pos[2]) * axb[2];
  C /= c_dot_axb;

  // the coefficients tell whether ppos lies within the cone defined by pos and
  // the adjacent edges

  if (A <= 0 && B <= 0 && C <= 0) {
    vec[0] = ppos[0] - pos[0];
    vec[1] = ppos[1] - pos[1];
    vec[2] = ppos[2] - pos[2];

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for cone at pos+a

  A = (ppos[0] - pos[0] - a[0]) * bxc[0] + (ppos[1] - pos[1] - a[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - a[0]) * axc[0] + (ppos[1] - pos[1] - a[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - a[0]) * axb[0] + (ppos[1] - pos[1] - a[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && B <= 0 && C <= 0) {
    vec[0] = ppos[0] - pos[0] - a[0];
    vec[1] = ppos[1] - pos[1] - a[1];
    vec[2] = ppos[2] - pos[2] - a[2];

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for cone at pos+b

  A = (ppos[0] - pos[0] - b[0]) * bxc[0] + (ppos[1] - pos[1] - b[1]) * bxc[1] +
      (ppos[2] - pos[2] - b[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - b[0]) * axc[0] + (ppos[1] - pos[1] - b[1]) * axc[1] +
      (ppos[2] - pos[2] - b[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - b[0]) * axb[0] + (ppos[1] - pos[1] - b[1]) * axb[1] +
      (ppos[2] - pos[2] - b[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && B >= 0 && C <= 0) {
    vec[0] = ppos[0] - pos[0] - b[0];
    vec[1] = ppos[1] - pos[1] - b[1];
    vec[2] = ppos[2] - pos[2] - b[2];

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for cone at pos+c

  A = (ppos[0] - pos[0] - c[0]) * bxc[0] + (ppos[1] - pos[1] - c[1]) * bxc[1] +
      (ppos[2] - pos[2] - c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - c[0]) * axc[0] + (ppos[1] - pos[1] - c[1]) * axc[1] +
      (ppos[2] - pos[2] - c[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - c[0]) * axb[0] + (ppos[1] - pos[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - c[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && B <= 0 && C >= 0) {
    vec[0] = ppos[0] - pos[0] - c[0];
    vec[1] = ppos[1] - pos[1] - c[1];
    vec[2] = ppos[2] - pos[2] - c[2];

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for cone at pos+a+b

  A = (ppos[0] - pos[0] - a[0] - b[0]) * bxc[0] +
      (ppos[1] - pos[1] - a[1] - b[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2] - b[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - a[0] - b[0]) * axc[0] +
      (ppos[1] - pos[1] - a[1] - b[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2] - b[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - a[0] - b[0]) * axb[0] +
      (ppos[1] - pos[1] - a[1] - b[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2] - b[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && B >= 0 && C <= 0) {
    vec[0] = ppos[0] - pos[0] - a[0] - b[0];
    vec[1] = ppos[1] - pos[1] - a[1] - b[1];
    vec[2] = ppos[2] - pos[2] - a[2] - b[2];

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for cone at pos+a+c

  A = (ppos[0] - pos[0] - a[0] - c[0]) * bxc[0] +
      (ppos[1] - pos[1] - a[1] - c[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2] - c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - a[0] - c[0]) * axc[0] +
      (ppos[1] - pos[1] - a[1] - c[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2] - c[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - a[0] - c[0]) * axb[0] +
      (ppos[1] - pos[1] - a[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2] - c[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && B <= 0 && C >= 0) {
    vec[0] = ppos[0] - pos[0] - a[0] - c[0];
    vec[1] = ppos[1] - pos[1] - a[1] - c[1];
    vec[2] = ppos[2] - pos[2] - a[2] - c[2];

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for cone at pos+a+c

  A = (ppos[0] - pos[0] - b[0] - c[0]) * bxc[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * bxc[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - b[0] - c[0]) * axc[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * axc[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - b[0] - c[0]) * axb[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && B >= 0 && C >= 0) {
    vec[0] = ppos[0] - pos[0] - b[0] - c[0];
    vec[1] = ppos[1] - pos[1] - b[1] - c[1];
    vec[2] = ppos[2] - pos[2] - b[2] - c[2];

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for cone at pos+a+b+c

  A = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * bxc[0] +
      (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * axc[0] +
      (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * axb[0] +
      (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && B >= 0 && C >= 0) {
    vec[0] = ppos[0] - pos[0] - a[0] - b[0] - c[0];
    vec[1] = ppos[1] - pos[1] - a[1] - b[1] - c[1];
    vec[2] = ppos[2] - pos[2] - a[2] - b[2] - c[2];

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos, a

  B = (ppos[0] - pos[0]) * axc[0] + (ppos[1] - pos[1]) * axc[1] +
      (ppos[2] - pos[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0]) * axb[0] + (ppos[1] - pos[1]) * axb[1] +
      (ppos[2] - pos[2]) * axb[2];
  C /= c_dot_axb;

  if (B <= 0 && C <= 0) {
    tmp = (ppos[0] - pos[0]) * a[0] + (ppos[1] - pos[1]) * a[1] +
          (ppos[2] - pos[2]) * a[2];
    tmp /= a[0] * a[0] + a[1] * a[1] + a[2] * a[2];

    vec[0] = ppos[0] - pos[0] - a[0] * tmp;
    vec[1] = ppos[1] - pos[1] - a[1] * tmp;
    vec[2] = ppos[2] - pos[2] - a[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos, b

  A = (ppos[0] - pos[0]) * bxc[0] + (ppos[1] - pos[1]) * bxc[1] +
      (ppos[2] - pos[2]) * bxc[2];
  A /= a_dot_bxc;
  C = (ppos[0] - pos[0]) * axb[0] + (ppos[1] - pos[1]) * axb[1] +
      (ppos[2] - pos[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && C <= 0) {
    tmp = (ppos[0] - pos[0]) * b[0] + (ppos[1] - pos[1]) * b[1] +
          (ppos[2] - pos[2]) * b[2];
    tmp /= b[0] * b[0] + b[1] * b[1] + b[2] * b[2];

    vec[0] = ppos[0] - pos[0] - b[0] * tmp;
    vec[1] = ppos[1] - pos[1] - b[1] * tmp;
    vec[2] = ppos[2] - pos[2] - b[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos, c

  A = (ppos[0] - pos[0]) * bxc[0] + (ppos[1] - pos[1]) * bxc[1] +
      (ppos[2] - pos[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0]) * axc[0] + (ppos[1] - pos[1]) * axc[1] +
      (ppos[2] - pos[2]) * axc[2];
  B /= b_dot_axc;

  if (A <= 0 && B <= 0) {
    tmp = (ppos[0] - pos[0]) * c[0] + (ppos[1] - pos[1]) * c[1] +
          (ppos[2] - pos[2]) * c[2];
    tmp /= c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

    vec[0] = ppos[0] - pos[0] - c[0] * tmp;
    vec[1] = ppos[1] - pos[1] - c[1] * tmp;
    vec[2] = ppos[2] - pos[2] - c[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+a, b

  A = (ppos[0] - pos[0] - a[0]) * bxc[0] + (ppos[1] - pos[1] - a[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2]) * bxc[2];
  A /= a_dot_bxc;
  C = (ppos[0] - pos[0] - a[0]) * axb[0] + (ppos[1] - pos[1] - a[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && C <= 0) {
    tmp = (ppos[0] - pos[0] - a[0]) * b[0] + (ppos[1] - pos[1] - a[1]) * b[1] +
          (ppos[2] - pos[2] - a[2]) * b[2];
    tmp /= b[0] * b[0] + b[1] * b[1] + b[2] * b[2];

    vec[0] = ppos[0] - pos[0] - a[0] - b[0] * tmp;
    vec[1] = ppos[1] - pos[1] - a[1] - b[1] * tmp;
    vec[2] = ppos[2] - pos[2] - a[2] - b[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+a, c

  A = (ppos[0] - pos[0] - a[0]) * bxc[0] + (ppos[1] - pos[1] - a[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - a[0]) * axc[0] + (ppos[1] - pos[1] - a[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2]) * axc[2];
  B /= b_dot_axc;

  if (A >= 0 && B <= 0) {
    tmp = (ppos[0] - pos[0] - a[0]) * c[0] + (ppos[1] - pos[1] - a[1]) * c[1] +
          (ppos[2] - pos[2] - a[2]) * c[2];
    tmp /= c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

    vec[0] = ppos[0] - pos[0] - a[0] - c[0] * tmp;
    vec[1] = ppos[1] - pos[1] - a[1] - c[1] * tmp;
    vec[2] = ppos[2] - pos[2] - a[2] - c[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+b+c, c

  A = (ppos[0] - pos[0] - b[0] - c[0]) * bxc[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * bxc[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - b[0] - c[0]) * axc[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * axc[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * axc[2];
  B /= b_dot_axc;

  if (A <= 0 && B >= 0) {
    tmp = (ppos[0] - pos[0] - b[0] - c[0]) * c[0] +
          (ppos[1] - pos[1] - b[1] - c[1]) * c[1] +
          (ppos[2] - pos[2] - b[2] - c[2]) * c[2];
    tmp /= c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

    vec[0] = ppos[0] - pos[0] - b[0] - c[0] - c[0] * tmp;
    vec[1] = ppos[1] - pos[1] - b[1] - c[1] - c[1] * tmp;
    vec[2] = ppos[2] - pos[2] - b[2] - c[2] - c[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+b+c, b

  A = (ppos[0] - pos[0] - b[0] - c[0]) * bxc[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * bxc[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * bxc[2];
  A /= a_dot_bxc;
  C = (ppos[0] - pos[0] - b[0] - c[0]) * axb[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * axb[2];
  C /= c_dot_axb;

  if (A <= 0 && C >= 0) {
    tmp = (ppos[0] - pos[0] - b[0] - c[0]) * b[0] +
          (ppos[1] - pos[1] - b[1] - c[1]) * b[1] +
          (ppos[2] - pos[2] - b[2] - c[2]) * b[2];
    tmp /= b[0] * b[0] + b[1] * b[1] + b[2] * b[2];

    vec[0] = ppos[0] - pos[0] - b[0] - c[0] - b[0] * tmp;
    vec[1] = ppos[1] - pos[1] - b[1] - c[1] - b[1] * tmp;
    vec[2] = ppos[2] - pos[2] - b[2] - c[2] - b[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+b+c, a

  B = (ppos[0] - pos[0] - b[0] - c[0]) * axc[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * axc[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - b[0] - c[0]) * axb[0] +
      (ppos[1] - pos[1] - b[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - b[2] - c[2]) * axb[2];
  C /= c_dot_axb;

  if (B >= 0 && C >= 0) {
    tmp = (ppos[0] - pos[0] - b[0] - c[0]) * a[0] +
          (ppos[1] - pos[1] - b[1] - c[1]) * a[1] +
          (ppos[2] - pos[2] - b[2] - c[2]) * a[2];
    tmp /= a[0] * a[0] + a[1] * a[1] + a[2] * a[2];

    vec[0] = ppos[0] - pos[0] - b[0] - c[0] - a[0] * tmp;
    vec[1] = ppos[1] - pos[1] - b[1] - c[1] - a[1] * tmp;
    vec[2] = ppos[2] - pos[2] - b[2] - c[2] - a[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+a+b, a

  B = (ppos[0] - pos[0] - a[0] - b[0]) * axc[0] +
      (ppos[1] - pos[1] - a[1] - b[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2] - b[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - a[0] - b[0]) * axb[0] +
      (ppos[1] - pos[1] - a[1] - b[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2] - b[2]) * axb[2];
  C /= c_dot_axb;

  if (B >= 0 && C <= 0) {
    tmp = (ppos[0] - pos[0] - a[0] - b[0]) * a[0] +
          (ppos[1] - pos[1] - a[1] - b[1]) * a[1] +
          (ppos[2] - pos[2] - a[2] - b[2]) * a[2];
    tmp /= a[0] * a[0] + a[1] * a[1] + a[2] * a[2];

    vec[0] = ppos[0] - pos[0] - a[0] - b[0] - a[0] * tmp;
    vec[1] = ppos[1] - pos[1] - a[1] - b[1] - a[1] * tmp;
    vec[2] = ppos[2] - pos[2] - a[2] - b[2] - a[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+a+b, c

  A = (ppos[0] - pos[0] - a[0] - b[0]) * bxc[0] +
      (ppos[1] - pos[1] - a[1] - b[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2] - b[2]) * bxc[2];
  A /= a_dot_bxc;
  B = (ppos[0] - pos[0] - a[0] - b[0]) * axc[0] +
      (ppos[1] - pos[1] - a[1] - b[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2] - b[2]) * axc[2];
  B /= b_dot_axc;

  if (A >= 0 && B >= 0) {
    tmp = (ppos[0] - pos[0] - a[0] - b[0]) * c[0] +
          (ppos[1] - pos[1] - a[1] - b[1]) * c[1] +
          (ppos[2] - pos[2] - a[2] - b[2]) * c[2];
    tmp /= c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

    vec[0] = ppos[0] - pos[0] - a[0] - b[0] - c[0] * tmp;
    vec[1] = ppos[1] - pos[1] - a[1] - b[1] - c[1] * tmp;
    vec[2] = ppos[2] - pos[2] - a[2] - b[2] - c[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+a+c, a

  B = (ppos[0] - pos[0] - a[0] - c[0]) * axc[0] +
      (ppos[1] - pos[1] - a[1] - c[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2] - c[2]) * axc[2];
  B /= b_dot_axc;
  C = (ppos[0] - pos[0] - a[0] - c[0]) * axb[0] +
      (ppos[1] - pos[1] - a[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2] - c[2]) * axb[2];
  C /= c_dot_axb;

  if (B <= 0 && C >= 0) {
    tmp = (ppos[0] - pos[0] - a[0] - c[0]) * a[0] +
          (ppos[1] - pos[1] - a[1] - c[1]) * a[1] +
          (ppos[2] - pos[2] - a[2] - c[2]) * a[2];
    tmp /= a[0] * a[0] + a[1] * a[1] + a[2] * a[2];

    vec[0] = ppos[0] - pos[0] - a[0] - c[0] - a[0] * tmp;
    vec[1] = ppos[1] - pos[1] - a[1] - c[1] - a[1] * tmp;
    vec[2] = ppos[2] - pos[2] - a[2] - c[2] - a[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for prism at edge pos+a+c, b

  A = (ppos[0] - pos[0] - a[0] - c[0]) * bxc[0] +
      (ppos[1] - pos[1] - a[1] - c[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2] - c[2]) * bxc[2];
  A /= a_dot_bxc;
  C = (ppos[0] - pos[0] - a[0] - c[0]) * axb[0] +
      (ppos[1] - pos[1] - a[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2] - c[2]) * axb[2];
  C /= c_dot_axb;

  if (A >= 0 && C >= 0) {
    tmp = (ppos[0] - pos[0] - a[0] - c[0]) * b[0] +
          (ppos[1] - pos[1] - a[1] - c[1]) * b[1] +
          (ppos[2] - pos[2] - a[2] - c[2]) * b[2];
    tmp /= b[0] * b[0] + b[1] * b[1] + b[2] * b[2];

    vec[0] = ppos[0] - pos[0] - a[0] - c[0] - b[0] * tmp;
    vec[1] = ppos[1] - pos[1] - a[1] - c[1] - b[1] * tmp;
    vec[2] = ppos[2] - pos[2] - a[2] - c[2] - b[2] * tmp;

    *dist =
        direction * sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    return 0;
  }

  // check for face with normal -axb

  *dist = (ppos[0] - pos[0]) * axb[0] + (ppos[1] - pos[1]) * axb[1] +
          (ppos[2] - pos[2]) * axb[2];
  *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
    *dist /= tmp;

    vec[0] = -*dist * axb[0] / tmp;
    vec[1] = -*dist * axb[1] / tmp;
    vec[2] = -*dist * axb[2] / tmp;

    *dist *= direction;

    return 0;
  }

  // calculate distance to face with normal axc

  *dist = (ppos[0] - pos[0]) * axc[0] + (ppos[1] - pos[1]) * axc[1] +
          (ppos[2] - pos[2]) * axc[2];

  if (*dist >= 0) {
    tmp = sqrt(axc[0] * axc[0] + axc[1] * axc[1] + axc[2] * axc[2]);
    *dist /= tmp;

    vec[0] = *dist * axc[0] / tmp;
    vec[1] = *dist * axc[1] / tmp;
    vec[2] = *dist * axc[2] / tmp;

    *dist *= direction;

    return 0;
  }

  // calculate distance to face with normal -bxc

  *dist = (ppos[0] - pos[0]) * bxc[0] + (ppos[1] - pos[1]) * bxc[1] +
          (ppos[2] - pos[2]) * bxc[2];
  *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(bxc[0] * bxc[0] + bxc[1] * bxc[1] + bxc[2] * bxc[2]);
    *dist /= tmp;

    vec[0] = -*dist * bxc[0] / tmp;
    vec[1] = -*dist * bxc[1] / tmp;
    vec[2] = -*dist * bxc[2] / tmp;

    *dist *= direction;

    return 0;
  }

  // calculate distance to face with normal axb

  *dist = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * axb[0] +
          (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * axb[1] +
          (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * axb[2];

  if (*dist >= 0) {
    tmp = sqrt(axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
    *dist /= tmp;

    vec[0] = *dist * axb[0] / tmp;
    vec[1] = *dist * axb[1] / tmp;
    vec[2] = *dist * axb[2] / tmp;

    *dist *= direction;

    return 0;
  }

  // calculate distance to face with normal -axc

  *dist = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * axc[0] +
          (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * axc[1] +
          (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * axc[2];
  *dist *= -1.;

  if (*dist >= 0) {
    tmp = sqrt(axc[0] * axc[0] + axc[1] * axc[1] + axc[2] * axc[2]);
    *dist /= tmp;

    vec[0] = -*dist * axc[0] / tmp;
    vec[1] = -*dist * axc[1] / tmp;
    vec[2] = -*dist * axc[2] / tmp;

    *dist *= direction;

    return 0;
  }

  // calculate distance to face with normal bxc

  *dist = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * bxc[0] +
          (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * bxc[1] +
          (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * bxc[2];

  if (*dist >= 0) {
    tmp = sqrt(bxc[0] * bxc[0] + bxc[1] * bxc[1] + bxc[2] * bxc[2]);
    *dist /= tmp;

    vec[0] = *dist * bxc[0] / tmp;
    vec[1] = *dist * bxc[1] / tmp;
    vec[2] = *dist * bxc[2] / tmp;

    *dist *= direction;

    return 0;
  }

  // ppos lies within rhomboid. Find nearest wall for interaction.

  // check for face with normal -axb

  *dist = (ppos[0] - pos[0]) * axb[0] + (ppos[1] - pos[1]) * axb[1] +
          (ppos[2] - pos[2]) * axb[2];
  *dist *= -1.;
  tmp = sqrt(axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
  *dist /= tmp;

  vec[0] = -*dist * axb[0] / tmp;
  vec[1] = -*dist * axb[1] / tmp;
  vec[2] = -*dist * axb[2] / tmp;

  *dist *= direction;

  // calculate distance to face with normal axc

  d = (ppos[0] - pos[0]) * axc[0] + (ppos[1] - pos[1]) * axc[1] +
      (ppos[2] - pos[2]) * axc[2];
  tmp = sqrt(axc[0] * axc[0] + axc[1] * axc[1] + axc[2] * axc[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = d * axc[0] / tmp;
    vec[1] = d * axc[1] / tmp;
    vec[2] = d * axc[2] / tmp;

    *dist = direction * d;
  }

  // calculate distance to face with normal -bxc

  d = (ppos[0] - pos[0]) * bxc[0] + (ppos[1] - pos[1]) * bxc[1] +
      (ppos[2] - pos[2]) * bxc[2];
  d *= -1.;
  tmp = sqrt(bxc[0] * bxc[0] + bxc[1] * bxc[1] + bxc[2] * bxc[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = -d * bxc[0] / tmp;
    vec[1] = -d * bxc[1] / tmp;
    vec[2] = -d * bxc[2] / tmp;

    *dist = direction * d;
  }

  // calculate distance to face with normal axb

  d = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * axb[0] +
      (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * axb[1] +
      (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * axb[2];
  tmp = sqrt(axb[0] * axb[0] + axb[1] * axb[1] + axb[2] * axb[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = d * axb[0] / tmp;
    vec[1] = d * axb[1] / tmp;
    vec[2] = d * axb[2] / tmp;

    *dist = direction * d;
  }

  // calculate distance to face with normal -axc

  d = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * axc[0] +
      (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * axc[1] +
      (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * axc[2];
  d *= -1.;
  tmp = sqrt(axc[0] * axc[0] + axc[1] * axc[1] + axc[2] * axc[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = -d * axc[0] / tmp;
    vec[1] = -d * axc[1] / tmp;
    vec[2] = -d * axc[2] / tmp;

    *dist = direction * d;
  }

  // calculate distance to face with normal bxc

  d = (ppos[0] - pos[0] - a[0] - b[0] - c[0]) * bxc[0] +
      (ppos[1] - pos[1] - a[1] - b[1] - c[1]) * bxc[1] +
      (ppos[2] - pos[2] - a[2] - b[2] - c[2]) * bxc[2];
  tmp = sqrt(bxc[0] * bxc[0] + bxc[1] * bxc[1] + bxc[2] * bxc[2]);
  d /= tmp;

  if (abs(d) < abs(*dist)) {
    vec[0] = d * bxc[0] / tmp;
    vec[1] = d * bxc[1] / tmp;
    vec[2] = d * bxc[2] / tmp;

    *dist = direction * d;
  }
  return 0;
}

}
