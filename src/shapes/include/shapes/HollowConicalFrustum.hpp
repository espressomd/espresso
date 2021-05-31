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
#include "utils/math/cylindrical_transformation_parameters.hpp"
#include <utils/Vector.hpp>

#include <list>
#include <memory>
#include <utility>

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
  HollowConicalFrustum(
      double const r1, double const r2, double const length,
      double const thickness, int const direction, double const central_angle,
      std::shared_ptr<Utils::CylindricalTransformationParameters>
          cyl_transform_params)
      : m_r1(r1), m_r2(r2), m_length(length), m_thickness(thickness),
        m_direction(direction), m_central_angle(central_angle),
        m_cyl_transform_params(std::move(cyl_transform_params)) {}

  void set_r1(double const radius) { m_r1 = radius; }
  void set_r2(double const radius) { m_r2 = radius; }
  void set_length(double const length) { m_length = length; }
  void set_thickness(double const thickness) { m_thickness = thickness; }
  void set_direction(int const dir) { m_direction = dir; }
  void set_central_angle(double const central_angle) {
    m_central_angle = central_angle;
  }

  /// Get radius 1 perpendicular to axis.
  double radius1() const { return m_r1; }
  /// Get radius 2 perpendicular to axis.
  double radius2() const { return m_r2; }
  /// Get length of the frustum (without thickness).
  double length() const { return m_length; }
  /// Get thickness of the frustum.
  double thickness() const { return m_thickness; }
  /// Get direction
  int direction() const { return m_direction; }
  /// Get central angle
  double central_angle() const { return m_central_angle; }

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
  double m_central_angle;
  std::shared_ptr<Utils::CylindricalTransformationParameters>
      m_cyl_transform_params;
};
} // namespace Shapes

#endif
