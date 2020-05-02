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
#ifndef CORE_BOX_GEOMETRY_HPP
#define CORE_BOX_GEOMETRY_HPP

#include <utils/Vector.hpp>

#include <bitset>
#include <cassert>

class BoxGeometry {
public:
  /** Flags for all three dimensions whether pbc are applied (default). */
  std::bitset<3> m_periodic = 0b111;
  /** Side lengths of the box */
  Utils::Vector3d m_length = {1, 1, 1};

  /**
   * @brief Set periodicity for direction
   *
   * @param coord The coordinate to set the periodicity for.
   * @param val True if this direction should be periodic.
   */
  void set_periodic(unsigned coord, bool val) { m_periodic.set(coord, val); }

  /**
   * @brief Check periodicity in direction.
   *
   * @param coord Direction to check
   * @return true iff periodic in direction.
   */
  constexpr bool periodic(unsigned coord) const {
    assert(coord <= 2);
    return m_periodic[coord];
  }

  /**
   * @brief Box length
   * @return Return vector of side-lengths of the box.
   */
  Utils::Vector3d const &length() const { return m_length; }

  /**
   * @brief Set box side lengths.
   * @param box_l Length that should be set.
   */
  void set_length(Utils::Vector3d const &box_l) { m_length = box_l; }

  /**
   * @brief Box volume
   * @return Return the volume of the box.
   */
  double volume() const { return m_length[0] * m_length[1] * m_length[2]; }
};

template <typename T> T get_mi_coord(T a, T b, T box_length, bool periodic) {
  auto const dx = a - b;

  if (periodic && (std::fabs(dx) > (0.5 * box_length))) {
    return dx - std::round(dx * (1. / box_length)) * box_length;
  }

  return dx;
}

template <typename T>
Utils::Vector<T, 3> get_mi_vector(const Utils::Vector<T, 3> &a,
                                  const Utils::Vector<T, 3> &b,
                                  const BoxGeometry &box) {
  return {get_mi_coord(a[0], b[0], box.length()[0], box.periodic(0)),
          get_mi_coord(a[1], b[1], box.length()[1], box.periodic(1)),
          get_mi_coord(a[2], b[2], box.length()[2], box.periodic(2))};
}

#endif // CORE_BOX_GEOMETRY_HPP
