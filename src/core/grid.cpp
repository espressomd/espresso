/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
/** \file
 *  Domain decomposition for parallel computing.
 *
 *  The corresponding header file is grid.hpp.
 */

#include "grid.hpp"

#include "communication.hpp"
#include "event.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <cmath>

BoxGeometry box_geo;
LocalBox local_geo;

LocalBox regular_decomposition(BoxGeometry const &box,
                               Utils::Vector3i const &node_index,
                               Utils::Vector3i const &node_grid) {

  auto const local_length = Utils::hadamard_division(box.length(), node_grid);
  auto const my_left = Utils::hadamard_product(node_index, local_length);

  Utils::Array<int, 6> boundaries;
  for (unsigned int dir = 0u; dir < 3u; dir++) {
    /* left boundary ? */
    boundaries[2u * dir] = (node_index[dir] == 0);
    /* right boundary ? */
    boundaries[2u * dir + 1u] = -(node_index[dir] + 1 == node_grid[dir]);
  }

  return {my_left, local_length, boundaries,
          CellStructureType::CELL_STRUCTURE_REGULAR};
}

void grid_changed_box_l(const BoxGeometry &box) {
  local_geo = regular_decomposition(box, ::communicator.calc_node_index(),
                                    ::communicator.node_grid);
}

void grid_changed_node_grid(bool update_box_geo) {
  ::communicator.init_comm_cart();

  if (update_box_geo) {
    grid_changed_box_l(box_geo);
  }
}

void set_node_grid(Utils::Vector3i const &value) {
  ::communicator.node_grid = value;
  on_node_grid_change();
}

void set_box_length(Utils::Vector3d const &value) {
  ::box_geo.set_length(value);
  on_boxl_change();
}
