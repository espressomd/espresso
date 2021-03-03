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

namespace Utils {

class CylTrafoParams {
public:
  explicit CylTrafoParams(
      Utils::Vector3d const &center = Utils::Vector3d{{0, 0, 0}},
      Utils::Vector3d const &axis = Utils::Vector3d{{0, 0, 1}},
      Utils::Vector3d const &orientation = Utils::Vector3d{{1, 0, 0}})
      : m_center(center), m_axis(axis), m_orientation(orientation) {
    check_valid();
  }

  Utils::Vector3d get_center() const { return m_center; }
  Utils::Vector3d get_axis() const { return m_axis; }
  Utils::Vector3d get_orientation() const { return m_orientation; }

  void set_center(Vector3d const &center) { m_center = center; }
  static void set_axis(Vector3d const &center) {
    throw std::runtime_error(
        "CylTrafoParams: Axis can only be set at construction.");
  }
  static void set_orientation(Vector3d const &center) {
    throw std::runtime_error(
        "CylTrafoParams: Orientation can only be set at construction.");
  }

private:
  void check_valid() {
    auto const eps = 10 * std::numeric_limits<double>::epsilon();
    if (Utils::abs(m_orientation * m_axis) > eps) {
      throw std::runtime_error("CylTrafoParams: Axis and orientation must be "
                               "orthogonal. Scalar product is " +
                               std::to_string(m_orientation * m_axis));
    }
    if (Utils::abs(m_axis.norm() - 1) > eps) {
      throw std::runtime_error(
          "CylTrafoParams: Axis must be normalized. Norm is " +
          std::to_string(m_axis.norm()));
    }
    if (Utils::abs(m_orientation.norm() - 1) > eps) {
      throw std::runtime_error(
          "CylTrafoParams: orientation must be normalized. Norm is " +
          std::to_string(m_orientation.norm()));
    }
  }

  Utils::Vector3d m_center;
  Utils::Vector3d m_axis;
  Utils::Vector3d m_orientation;
};

} // namespace Utils

#endif // ESPRESSO_CYLINDER_TRANSFORMATION_PARAMETERS_HPP