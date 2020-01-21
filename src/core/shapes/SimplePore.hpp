/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#ifndef SHAPES_SIMPLE_PORE_HPP
#define SHAPES_SIMPLE_PORE_HPP

#include "Shape.hpp"
#include <utils/Vector.hpp>

namespace Shapes {
class SimplePore : public Shape {
  double m_rad;
  double m_length;
  double m_smoothing_rad;
  Utils::Vector3d m_center;
  Utils::Vector3d m_axis;

  /* Center of smoothing circle */
  double c_r;
  double c_z;
  double m_half_length;
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

    c_r = m_rad + m_smoothing_rad;
    c_z = m_half_length - m_smoothing_rad;
  }

  std::pair<double, double> dist_half_pore(double r, double z) const;

public:
  SimplePore() : m_axis({1., 0., 0}) { precalc(); }

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

  double smoothing_radius() const { return m_smoothing_rad; }
  void set_smoothing_radius(double const &smoothing_radius) {
    m_smoothing_rad = smoothing_radius;
    precalc();
  }

  Utils::Vector3d const &axis() const { return m_axis; }
  void set_axis(Utils::Vector3d const &axis) {
    m_axis = axis;
    precalc();
  }

  Utils::Vector3d &center() { return m_center; }

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;
};
} // namespace Shapes

#endif
