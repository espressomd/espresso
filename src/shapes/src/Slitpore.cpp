/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
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

#include <shapes/Slitpore.hpp>

#include <utils/math/sqr.hpp>

#include <cmath>

using namespace std;

namespace Shapes {
void Slitpore::calculate_dist(const Utils::Vector3d &pos, double &dist,
                              Utils::Vector3d &vec) const {
  // the left circles
  Utils::Vector2d c11 = {dividing_plane() - m_pore_width / 2 -
                             m_upper_smoothing_radius,
                         m_pore_mouth - m_upper_smoothing_radius};
  Utils::Vector2d c12 = {
      dividing_plane() - m_pore_width / 2 + m_lower_smoothing_radius,
      m_pore_mouth - m_pore_length + m_lower_smoothing_radius};
  // the right circles
  Utils::Vector2d c21 = {dividing_plane() + m_pore_width / 2 +
                             m_upper_smoothing_radius,
                         m_pore_mouth - m_upper_smoothing_radius};
  Utils::Vector2d c22 = {
      dividing_plane() + m_pore_width / 2 - m_lower_smoothing_radius,
      m_pore_mouth - m_pore_length + m_lower_smoothing_radius};

  if (pos[2] > m_pore_mouth + m_channel_width / 2) {
    // Feel the upper wall
    dist = m_pore_mouth + m_channel_width - pos[2];
    vec[0] = vec[1] = 0;
    vec[2] = -dist;
    return;
  }

  if (pos[0] < c11[0] || pos[0] > c21[0]) {
    // Feel the lower wall of the channel
    dist = pos[2] - m_pore_mouth;
    vec[0] = vec[1] = 0;
    vec[2] = dist;
    return;
  }

  if (pos[2] > c11[1]) {
    // Feel the upper smoothing
    if (pos[0] < dividing_plane()) {
      dist = sqrt(Utils::sqr(c11[0] - pos[0]) + Utils::sqr(c11[1] - pos[2])) -
             m_upper_smoothing_radius;
      vec[0] = -(c11[0] - pos[0]) * (dist) / (dist + m_upper_smoothing_radius);
      vec[1] = 0;
      vec[2] = -(c11[1] - pos[2]) * (dist) / (dist + m_upper_smoothing_radius);
      return;
    }
    dist = sqrt(Utils::sqr(c21[0] - pos[0]) + Utils::sqr(c21[1] - pos[2])) -
           m_upper_smoothing_radius;
    vec[0] = -(c21[0] - pos[0]) * (dist) / (dist + m_upper_smoothing_radius);
    vec[1] = 0;
    vec[2] = -(c21[1] - pos[2]) * (dist) / (dist + m_upper_smoothing_radius);
    return;
  }

  if (pos[2] > c12[1]) {
    // Feel the pore wall
    if (pos[0] < dividing_plane()) {
      dist = pos[0] - (dividing_plane() - m_pore_width / 2);
      vec[0] = dist;
      vec[1] = vec[2] = 0;
      return;
    }
    dist = (dividing_plane() + m_pore_width / 2) - pos[0];
    vec[0] = -dist;
    vec[1] = vec[2] = 0;
    return;
  }

  if (pos[0] > c12[0] && pos[0] < c22[0]) {
    // Feel the pore end wall
    dist = pos[2] - (m_pore_mouth - m_pore_length);
    vec[0] = vec[1] = 0;
    vec[2] = dist;
    return;
  }

  // Feel the lower smoothing
  if (pos[0] < dividing_plane()) {
    dist = -sqrt(Utils::sqr(c12[0] - pos[0]) + Utils::sqr(c12[1] - pos[2])) +
           m_lower_smoothing_radius;
    vec[0] = (c12[0] - pos[0]) * (dist) / (-dist + m_lower_smoothing_radius);
    vec[1] = 0;
    vec[2] = (c12[1] - pos[2]) * (dist) / (-dist + m_lower_smoothing_radius);
    return;
  }
  dist = -sqrt(Utils::sqr(c22[0] - pos[0]) + Utils::sqr(c22[1] - pos[2])) +
         m_lower_smoothing_radius;
  vec[0] = (c22[0] - pos[0]) * (dist) / (-dist + m_lower_smoothing_radius);
  vec[1] = 0;
  vec[2] = (c22[1] - pos[2]) * (dist) / (-dist + m_lower_smoothing_radius);
}
} // namespace Shapes
