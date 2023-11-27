/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "cell_system/CellStructureType.hpp"

#include <utils/Array.hpp>
#include <utils/Vector.hpp>

class LocalBox {
  Utils::Vector3d m_local_box_l = {1., 1., 1.};
  Utils::Vector3d m_lower_corner = {0., 0., 0.};
  Utils::Vector3d m_upper_corner = {1., 1., 1.};
  Utils::Array<int, 6> m_boundaries = {};
  CellStructureType m_cell_structure_type;

public:
  LocalBox() = default;
  LocalBox(Utils::Vector3d const &lower_corner,
           Utils::Vector3d const &local_box_length,
           Utils::Array<int, 6> const &boundaries,
           CellStructureType const cell_structure_type)
      : m_local_box_l(local_box_length), m_lower_corner(lower_corner),
        m_upper_corner(lower_corner + local_box_length),
        m_boundaries(boundaries), m_cell_structure_type(cell_structure_type) {}

  /** Left (bottom, front) corner of this nodes local box. */
  auto const &my_left() const { return m_lower_corner; }
  /** Right (top, back) corner of this nodes local box. */
  auto const &my_right() const { return m_upper_corner; }
  /** Dimensions of the box a single node is responsible for. */
  auto const &length() const { return m_local_box_l; }
  /** @brief Boundary information for the local box.
   *
   * This returns for each of the faces of the local box if
   * it is a boundary of the simulation box. The format is
   * as follows:
   *  (x low, x high, y low, y high, z low, z high).
   *
   * @return Array with boundary information.
   */
  auto const &boundary() const { return m_boundaries; }

  /** Return cell structure type. */
  auto const &cell_structure_type() const { return m_cell_structure_type; }

  /** Set cell structure type. */
  void set_cell_structure_type(CellStructureType cell_structure_type) {
    m_cell_structure_type = cell_structure_type;
  }

  static LocalBox make_regular_decomposition(Utils::Vector3d const &box_l,
                                             Utils::Vector3i const &node_index,
                                             Utils::Vector3i const &node_grid) {

    auto const local_length = Utils::hadamard_division(box_l, node_grid);
    auto const my_left = Utils::hadamard_product(node_index, local_length);

    decltype(LocalBox::m_boundaries) boundaries;
    for (unsigned int dir = 0u; dir < 3u; dir++) {
      /* left boundary ? */
      boundaries[2u * dir] = (node_index[dir] == 0);
      /* right boundary ? */
      boundaries[2u * dir + 1u] = -(node_index[dir] + 1 == node_grid[dir]);
    }

    return {my_left, local_length, boundaries, CellStructureType::REGULAR};
  }
};
