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

#ifndef __SPHEROCYLINDER_HPP
#define __SPHEROCYLINDER_HPP

#include "Shape.hpp"
#include <utils/Vector.hpp>

namespace Shapes {
class SpheroCylinder : public Shape {
public:
  /** center of the cylinder. */
  Utils::Vector3d m_center;
  /** Axis of the cylinder. */
  Utils::Vector3d m_axis;
  /** cylinder radius. */
  double m_rad;
  /** cylinder length. */
  double m_length;

  /** Center of smoothing circle */
  double m_half_length;
  /** direction -1: inside, +1 outside */
  double m_direction;
  /** Unit vector in z direction */
  Utils::Vector3d e_z;

  /** Alternative e_r for corner case */
  Utils::Vector3d e_r_axis;

  /** @brief Calculate derived parameters. */
  void precalc() {
    m_half_length = 0.5 * m_length;

    e_z = m_axis / m_axis.norm();

    /* Find a vector orthogonal to e_z, since {1,0,0} and
       {0,1,0} are independent, e_z can not be parallel to both
       of them. Then we can do Gram-Schmidt */
    if ((Utils::Vector3d{1., 0., 0} * e_z) < 1.)
      e_r_axis =
          Utils::Vector3d{1., 0., 0} - (e_z * Utils::Vector3d{1., 0., 0}) * e_z;
    else
      e_r_axis =
          Utils::Vector3d{0., 1., 0} - (e_z * Utils::Vector3d{0., 1., 0}) * e_z;

    e_r_axis.normalize();
  }

public:
  SpheroCylinder()
      : m_center({0.0, 0.0, 0.0}), m_axis({1.0, 0.0, 0.0}), m_rad(0),
        m_length(0.0) {
    precalc();
  }

  double radius() const { return m_rad; }
  void set_radius(double const &radius) {
    m_rad = radius;
    precalc();
  }

  double length() const { return m_length; }
  void set_length(double const &length) {
    m_length = length;
    precalc();
  }

  Utils::Vector3d const &axis() const { return m_axis; }
  void set_axis(Utils::Vector3d const &axis) {
    m_axis = axis;
    precalc();
  }

  Utils::Vector3d &center() { return m_center; }
  double &direction() { return m_direction; }

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;
};
} // namespace Shapes

#endif
