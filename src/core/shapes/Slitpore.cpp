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

#include "Slitpore.hpp"
#include "utils.hpp"
/* For box_l */
#include "grid.hpp"

#include <cmath>

using namespace std;

namespace Shapes {
int Slitpore::calculate_dist(const double *ppos, double *dist,
                             double *vec) const {
  // the left circles
  double box_l_x = box_l[0];
  double c11[2] = {box_l_x / 2 - m_pore_width / 2 - m_upper_smoothing_radius,
                   m_pore_mouth - m_upper_smoothing_radius};
  double c12[2] = {box_l_x / 2 - m_pore_width / 2 + m_lower_smoothing_radius,
                   m_pore_mouth - m_pore_length + m_lower_smoothing_radius};
  // the right circles
  double c21[2] = {box_l_x / 2 + m_pore_width / 2 + m_upper_smoothing_radius,
                   m_pore_mouth - m_upper_smoothing_radius};
  double c22[2] = {box_l_x / 2 + m_pore_width / 2 - m_lower_smoothing_radius,
                   m_pore_mouth - m_pore_length + m_lower_smoothing_radius};

  if (ppos[2] > m_pore_mouth + m_channel_width / 2) {
    //    printf("upper wall\n");
    // Feel the upper wall
    *dist = m_pore_mouth + m_channel_width - ppos[2];
    vec[0] = vec[1] = 0;
    vec[2] = -*dist;
    return 0;
  }

  if (ppos[0] < c11[0] || ppos[0] > c21[0]) {
    // Feel the lower wall of the channel
    *dist = ppos[2] - m_pore_mouth;
    vec[0] = vec[1] = 0;
    vec[2] = *dist;
    return 0;
  }

  if (ppos[2] > c11[1]) {
    // Feel the upper smoothing
    if (ppos[0] < box_l_x / 2) {
      *dist =
          sqrt(Utils::sqr(c11[0] - ppos[0]) + Utils::sqr(c11[1] - ppos[2])) -
          m_upper_smoothing_radius;
      vec[0] =
          -(c11[0] - ppos[0]) * (*dist) / (*dist + m_upper_smoothing_radius);
      vec[1] = 0;
      vec[2] =
          -(c11[1] - ppos[2]) * (*dist) / (*dist + m_upper_smoothing_radius);
      return 0;
    } else {
      *dist =
          sqrt(Utils::sqr(c21[0] - ppos[0]) + Utils::sqr(c21[1] - ppos[2])) -
          m_upper_smoothing_radius;
      vec[0] =
          -(c21[0] - ppos[0]) * (*dist) / (*dist + m_upper_smoothing_radius);
      vec[1] = 0;
      vec[2] =
          -(c21[1] - ppos[2]) * (*dist) / (*dist + m_upper_smoothing_radius);
      return 0;
    }
  }

  if (ppos[2] > c12[1]) {
    // Feel the pore wall
    if (ppos[0] < box_l_x / 2) {
      *dist = ppos[0] - (box_l_x / 2 - m_pore_width / 2);
      vec[0] = *dist;
      vec[1] = vec[2] = 0;
      return 0;
    } else {
      *dist = (box_l_x / 2 + m_pore_width / 2) - ppos[0];
      vec[0] = -*dist;
      vec[1] = vec[2] = 0;
      return 0;
    }
  }

  if (ppos[0] > c12[0] && ppos[0] < c22[0]) {
    //    printf("pore end\n");
    // Feel the pore end wall
    *dist = ppos[2] - (m_pore_mouth - m_pore_length);
    vec[0] = vec[1] = 0;
    vec[2] = *dist;
    return 0;
  }
  // Else
  // Feel the lower smoothing
  if (ppos[0] < box_l_x / 2) {
    *dist = -sqrt(Utils::sqr(c12[0] - ppos[0]) + Utils::sqr(c12[1] - ppos[2])) +
            m_lower_smoothing_radius;
    vec[0] = (c12[0] - ppos[0]) * (*dist) / (-*dist + m_lower_smoothing_radius);
    vec[1] = 0;
    vec[2] = (c12[1] - ppos[2]) * (*dist) / (-*dist + m_lower_smoothing_radius);
    return 0;
  } else {
    *dist = -sqrt(Utils::sqr(c22[0] - ppos[0]) + Utils::sqr(c22[1] - ppos[2])) +
            m_lower_smoothing_radius;
    vec[0] = (c22[0] - ppos[0]) * (*dist) / (-*dist + m_lower_smoothing_radius);
    vec[1] = 0;
    vec[2] = (c22[1] - ppos[2]) * (*dist) / (-*dist + m_lower_smoothing_radius);
    return 0;
  }

  return 0;
}
} // namespace Shapes
