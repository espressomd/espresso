/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef ESPRESSO_CYLINDER_TRANSFORMATION_PARAMETERS_HPP
#define ESPRESSO_CYLINDER_TRANSFORMATION_PARAMETERS_HPP

#include <stdexcept>
#include <string>

#include <utils/math/abs.hpp>
#include <utils/math/orthonormal_vec.hpp>

namespace Utils {

/**
 * @brief A class to hold and validate parameters for a cylindrical coordinate
 * transformations.
 *
 * @param center The origin of the cylindrical coordinates.
 * @param axis The "z"-axis. Must be normalized.
 * @param orientation The axis along which phi = 0. Must be normalized and
 * orthogonal to axis.
 */
class CylindricalTransformationParameters {
public:
  CylindricalTransformationParameters() = default;
  CylindricalTransformationParameters(Utils::Vector3d const &center,
                                      Utils::Vector3d const &axis,
                                      Utils::Vector3d const &orientation)
      : m_center(center), m_axis(axis), m_orientation(orientation) {
    validate();
  }
  /**
   * @brief if you only provide center and axis, an orientation will be
   * generated automatically such that it is orthogonal to axis
   */
  CylindricalTransformationParameters(Utils::Vector3d const &center,
                                      Utils::Vector3d const &axis)
      : m_center(center), m_axis(axis),
        m_orientation(Utils::calc_orthonormal_vector(axis)) {}

  Utils::Vector3d center() const { return m_center; }
  Utils::Vector3d axis() const { return m_axis; }
  Utils::Vector3d orientation() const { return m_orientation; }

private:
  void validate() const {
    auto constexpr eps = 10 * std::numeric_limits<double>::epsilon();
    if (Utils::abs(m_orientation * m_axis) > eps) {
      throw std::runtime_error(
          "CylindricalTransformationParameters: Axis and orientation must be "
          "orthogonal. Scalar product is " +
          std::to_string(m_orientation * m_axis));
    }
    if (Utils::abs(m_axis.norm() - 1) > eps) {
      throw std::runtime_error("CylindricalTransformationParameters: Axis must "
                               "be normalized. Norm is " +
                               std::to_string(m_axis.norm()));
    }
    if (Utils::abs(m_orientation.norm() - 1) > eps) {
      throw std::runtime_error("CylindricalTransformationParameters: "
                               "orientation must be normalized. Norm is " +
                               std::to_string(m_orientation.norm()));
    }
  }

  const Utils::Vector3d m_center{};
  const Utils::Vector3d m_axis{0, 0, 1};
  const Utils::Vector3d m_orientation{1, 0, 0};
};

} // namespace Utils

#endif // ESPRESSO_CYLINDER_TRANSFORMATION_PARAMETERS_HPP
