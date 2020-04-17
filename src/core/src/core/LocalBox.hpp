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
#ifndef ESPRESSO_LOCALBOX_HPP
#define ESPRESSO_LOCALBOX_HPP

#include <utils/Vector.hpp>

template <class T> class LocalBox {
  Utils::Vector<T, 3> m_local_box_l = {1, 1, 1};
  Utils::Vector<T, 3> m_lower_corner = {0, 0, 0};
  Utils::Vector<T, 3> m_upper_corner = {1, 1, 1};
  Utils::Array<int, 6> m_boundaries = {};

public:
  LocalBox() = default;
  LocalBox(Utils::Vector<T, 3> const &lower_corner,
           Utils::Vector<T, 3> const &local_box_length,
           Utils::Array<int, 6> const &boundaries)
      : m_local_box_l(local_box_length), m_lower_corner(lower_corner),
        m_upper_corner(lower_corner + local_box_length),
        m_boundaries(boundaries) {}

  /** Left (bottom, front) corner of this nodes local box. */
  Utils::Vector<T, 3> const &my_left() const { return m_lower_corner; }
  /** Right (top, back) corner of this nodes local box. */
  Utils::Vector<T, 3> const &my_right() const { return m_upper_corner; }
  /** Dimensions of the box a single node is responsible for. */
  Utils::Vector<T, 3> const &length() const { return m_local_box_l; }
  /** @brief Boundary information for the local box.
   *
   * This returns for each of the faces of the local box if
   * it is a boundary of the simulation box. The format is
   * as follows:
   *  (x low, x high, y low, y high, z low, z high).
   *
   * @return Array with boundary information.
   */
  Utils::Array<int, 6> const &boundary() const { return m_boundaries; }
};

#endif
