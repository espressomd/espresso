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

#ifndef SHAPES_CONICAL_FRUSTUM_HPP
#define SHAPES_CONICAL_FRUSTUM_HPP

#include "Shape.hpp"
#include <utils/Vector.hpp>

namespace Shapes {

/**
 * <pre>
 * -------->r1
 * ^ a      \
 * | x       \
 * | i        \
 * | s         \                       *pos
 * * center     \
 * |             \
 * |              \
 * |               \
 * |                \
 * ----------------->r2
 * </pre>
 * @brief Conical frustum shape with rounded corners and finite thickness.
 *
 * \image html conical_frustum.png width=800px
 */
class HollowConicalFrustum : public Shape {
public:
  HollowConicalFrustum()
      : m_r1(0.0), m_r2(0.0), m_length(0.0), m_thickness(0.0),
        m_direction(1), m_center{Utils::Vector3d{}}, m_axis{Utils::Vector3d{
                                                         0, 0, 1}} {}

  void set_r1(double radius) { m_r1 = radius; }
  void set_r2(double radius) { m_r2 = radius; }
  void set_length(double length) { m_length = length; }
  void set_thickness(double thickness) { m_thickness = thickness; }
  void set_direction(int dir) { m_direction = dir; }
  void set_axis(Utils::Vector3d const &axis) { m_axis = axis; }

  void set_center(Utils::Vector3d const &center) { m_center = center; }

  /// Get radius 1 perpendicular to axis.
  double radius1() const { return m_r1; }
  /// Get radius 2 perpendicular to axis.
  double radius2() const { return m_r2; }
  /// Get length of the frustum (modulo thickness).
  double length() const { return m_length; }
  /// Get thickness of the frustum.
  double thickness() const { return m_thickness; }
  /// Get the direction of the shape. If -1, distance is positive within the
  /// enclosed volume of the frustum.
  int direction() const { return m_direction; }
  /// Get center of symmetry.
  Utils::Vector3d const &center() const { return m_center; }
  /// Get symmetry axis.
  Utils::Vector3d const &axis() const { return m_axis; }
  /**
   * @brief Calculate the distance vector and its norm between a given position
   * and the cone.
   * @param[in] pos Position for which distance to the cone is calculated.
   * @param[out] dist Distance between cone and \p pos.
   * @param[out] vec Distance vector (\p dist = || \p vec ||).
   */
  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;

private:
  double m_r1;
  double m_r2;
  double m_length;
  double m_thickness;
  int m_direction;
  Utils::Vector3d m_center;
  Utils::Vector3d m_axis;
};
} // namespace Shapes

#endif
