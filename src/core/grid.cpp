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
#include <utils/mpi/cart_comm.hpp>

#include <mpi.h>

#include <algorithm>
#include <cmath>

BoxGeometry box_geo;
LocalBox<double> local_geo;

Utils::Vector3i node_grid{};

void init_node_grid() { grid_changed_n_nodes(); }

Utils::Vector3i calc_node_pos(const boost::mpi::communicator &comm) {
  return Utils::Mpi::cart_coords<3>(comm, comm.rank());
}

Utils::Vector<int, 6>
calc_node_neighbors(const boost::mpi::communicator &comm) {
  return Utils::Mpi::cart_neighbors<3>(comm);
}

LocalBox<double> regular_decomposition(const BoxGeometry &box,
                                       Utils::Vector3i const &node_pos,
                                       Utils::Vector3i const &node_grid_par) {
  Utils::Vector3d local_length;
  Utils::Vector3d my_left;

  for (unsigned int i = 0; i < 3; i++) {
    local_length[i] = box.length()[i] / node_grid_par[i];
    my_left[i] = node_pos[i] * local_length[i];
  }

  Utils::Array<int, 6> boundaries;
  for (unsigned int dir = 0; dir < 3; dir++) {
    /* left boundary ? */
    boundaries[2 * dir] = (node_pos[dir] == 0);
    /* right boundary ? */
    boundaries[2 * dir + 1] = -(node_pos[dir] == node_grid_par[dir] - 1);
  }

  return {my_left, local_length, boundaries,
          CellStructureType::CELL_STRUCTURE_REGULAR};
}

void grid_changed_box_l(const BoxGeometry &box) {
  local_geo = regular_decomposition(box, calc_node_pos(comm_cart), node_grid);
}

void grid_changed_n_nodes() {
  comm_cart =
      Utils::Mpi::cart_create(comm_cart, node_grid, /* reorder */ false);

  this_node = comm_cart.rank();

  calc_node_neighbors(comm_cart);

  grid_changed_box_l(box_geo);
}

void set_node_grid(Utils::Vector3i const &value) {
  ::node_grid = value;
  on_node_grid_change();
}

void set_box_length(Utils::Vector3d const &value) {
  ::box_geo.set_length(value);
  on_boxl_change();
}
