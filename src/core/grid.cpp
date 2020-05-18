/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "cells.hpp"
#include "communication.hpp"
#include "global.hpp"

#include <boost/algorithm/clamp.hpp>
#include <mpi.h>
#include <utils/mpi/cart_comm.hpp>

/**********************************************
 * variables
 **********************************************/

BoxGeometry box_geo;
LocalBox<double> local_geo;

Utils::Vector3i node_grid{};

/************************************************************/

void init_node_grid() { grid_changed_n_nodes(); }

int map_position_node_array(const Utils::Vector3d &pos) {
  auto const f_pos = folded_position(pos, box_geo);

  Utils::Vector3i im;
  for (int i = 0; i < 3; i++) {
    im[i] = std::floor(f_pos[i] / local_geo.length()[i]);
    im[i] = boost::algorithm::clamp(im[i], 0, node_grid[i] - 1);
  }

  return Utils::Mpi::cart_rank(comm_cart, im);
}

Utils::Vector3i calc_node_pos(const boost::mpi::communicator &comm) {
  return Utils::Mpi::cart_coords<3>(comm, comm.rank());
}

Utils::Vector<int, 6>
calc_node_neighbors(const boost::mpi::communicator &comm) {
  return Utils::Mpi::cart_neighbors<3>(comm);
}

LocalBox<double> regular_decomposition(const BoxGeometry &box,
                                       Utils::Vector3i const &node_pos,
                                       Utils::Vector3i const &node_grid) {
  Utils::Vector3d local_length;
  Utils::Vector3d my_left;

  for (int i = 0; i < 3; i++) {
    local_length[i] = box.length()[i] / node_grid[i];
    my_left[i] = node_pos[i] * local_length[i];
  }

  Utils::Array<int, 6> boundaries;
  for (int dir = 0; dir < 3; dir++) {
    /* left boundary ? */
    boundaries[2 * dir] = (node_pos[dir] == 0);
    /* right boundary ? */
    boundaries[2 * dir + 1] = -(node_pos[dir] == node_grid[dir] - 1);
  }

  return {my_left, local_length, boundaries};
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

void rescale_boxl(int dir, double d_new) {
  double scale =
      (dir - 3) ? d_new / box_geo.length()[dir] : d_new / box_geo.length()[0];

  /* If shrinking, rescale the particles first. */
  if (scale <= 1.) {
    mpi_rescale_particles(dir, scale);
  }

  if (dir < 3) {
    auto box_l = box_geo.length();
    box_l[dir] = d_new;
    box_geo.set_length(box_l);
  } else {
    box_geo.set_length({d_new, d_new, d_new});
  }

  mpi_bcast_parameter(FIELD_BOXL);

  if (scale > 1.) {
    mpi_rescale_particles(dir, scale);
  }
}
