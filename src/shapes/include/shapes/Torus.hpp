/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include "Shape.hpp"

#include <utils/Vector.hpp>

namespace Shapes {
class Torus : public Shape {
public:
  /** center of the shape. */
  Utils::Vector3d m_center;
  /** Normal axis of the shape. */
  Utils::Vector3d m_normal;
  /** radius. */
  double m_rad;
  /** tube radius. */
  double m_tube_rad;
  /** direction -1: inside, +1 outside */
  double m_direction;

  /** Unit vector in z direction */
  Utils::Vector3d e_z;

  /** Alternative e_r for the on-axis corner case */
  Utils::Vector3d e_r_axis;

  /** @brief Calculate derived parameters. */
  void precalc() {
    e_z = m_normal / m_normal.norm();

    /* Find a vector orthogonal to e_z via Gram-Schmidt.
       Since {1,0,0} and {0,1,0} are independent, e_z cannot be
       parallel to both of them. */
    if (e_z[0] < 1.) {
      Utils::Vector3d constexpr u_x{{1., 0., 0.}};
      e_r_axis = u_x - (e_z * u_x) * e_z;
    } else {
      Utils::Vector3d constexpr u_y{{0., 1., 0.}};
      e_r_axis = u_y - (e_z * u_y) * e_z;
    }
    e_r_axis.normalize();
  }

public:
  Torus()
      : m_center({0.0, 0.0, 0.0}), m_normal({1.0, 0.0, 0.0}), m_rad(0.0),
        m_tube_rad(0.0), m_direction(1.0) {
    precalc();
  }

  double radius() const { return m_rad; }
  void set_radius(double const &radius) {
    m_rad = radius;
    precalc();
  }

  double tube_radius() const { return m_tube_rad; }
  void set_tube_radius(double const &tube_rad) {
    m_tube_rad = tube_rad;
    precalc();
  }

  Utils::Vector3d const &normal() const { return m_normal; }
  void set_normal(Utils::Vector3d const &normal) {
    m_normal = normal;
    precalc();
  }

  Utils::Vector3d &center() { return m_center; }
  double &direction() { return m_direction; }

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;
};
} // namespace Shapes
