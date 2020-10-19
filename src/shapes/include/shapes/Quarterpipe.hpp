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

#ifndef SHAPES_QUARTERPIPE_HPP
#define SHAPES_QUARTERPIPE_HPP

#include "Shape.hpp"
#include <utils/Vector.hpp>
#include <utils/constants.hpp>

namespace Shapes {

/**
 * @brief Quarterpipe shape based on intersecting a cuboid
 * with a cylinder
 *
 * \image html quarterpipe_detailed.png width=1000px
 */

class Quarterpipe : public Shape {
public:
  Quarterpipe()
      : m_center{Utils::Vector3d{}}, m_axis{Utils::Vector3d{}},
        m_orientation{Utils::Vector3d{}}, m_radius(0.), m_height(0.),
        m_delta_phi(0.5 * Utils::pi()) {}

  void set_center(Utils::Vector3d const &center) { m_center = center; }
  void set_axis(Utils::Vector3d const &axis) {
    auto axis_tmp(axis);
    m_axis = axis_tmp.normalize();
  }
  void set_orientation(Utils::Vector3d const &orientation) {
    auto orientation_tmp(orientation);
    m_orientation = orientation_tmp.normalize();
  }
  void set_radius(double const radius) { m_radius = radius; }
  void set_height(double const height) { m_height = height; }

  Utils::Vector3d const &center() const { return m_center; }
  Utils::Vector3d const &axis() const { return m_axis; }
  Utils::Vector3d const &orientation() const { return m_orientation; }
  double radius() const { return m_radius; }
  double height() const { return m_height; }

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;

private:
  Utils::Vector3d m_center;
  Utils::Vector3d m_axis;
  Utils::Vector3d m_orientation;
  double m_radius;
  double m_height;
  double m_delta_phi;
};
} // namespace Shapes

#endif
