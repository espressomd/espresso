/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *  Domain decomposition for parallel computing.
 *
 *  The corresponding header file is grid.hpp.
 */

#include "grid.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "debug.hpp"
#include "global.hpp"

#include <boost/algorithm/clamp.hpp>
#include <mpi.h>

#include <boost/range/algorithm/min_element.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/**********************************************
 * variables
 **********************************************/

BoxGeometry box_geo;
LocalBox<double> local_geo;

Utils::Vector3i node_grid{};
Utils::Vector3i node_pos = {-1, -1, -1};
Utils::Vector<int, 6> node_neighbors{};

double min_box_l;
double min_local_box_l;

/************************************************************/

void init_node_grid() {
  grid_changed_n_nodes();
  cells_on_geometry_change(CELL_FLAG_GRIDCHANGED);
}

int map_position_node_array(const Utils::Vector3d &pos) {
  auto const f_pos = folded_position(pos, box_geo);

  Utils::Vector3i im;
  for (int i = 0; i < 3; i++) {
    im[i] = std::floor(f_pos[i] / local_geo.length()[i]);
    im[i] = boost::algorithm::clamp(im[i], 0, node_grid[i] - 1);
  }

  auto const node = map_array_node(im);
  return node;
}

int calc_node_neighbors(int node) {

  int dir, neighbor_count;

  map_node_array(node, node_pos.data());
  for (dir = 0; dir < 3; dir++) {
    int buf;

    MPI_Cart_shift(comm_cart, dir, -1, &buf, &(node_neighbors[2 * dir]));
    MPI_Cart_shift(comm_cart, dir, 1, &buf, &(node_neighbors[2 * dir + 1]));

    /* left boundary ? */
    if (node_pos[dir] == 0) {
      local_geo.boundary_[2 * dir] = 1;
    } else {
      local_geo.boundary_[2 * dir] = 0;
    }
    /* right boundary ? */
    if (node_pos[dir] == node_grid[dir] - 1) {
      local_geo.boundary_[2 * dir + 1] = -1;
    } else {
      local_geo.boundary_[2 * dir + 1] = 0;
    }
  }

  neighbor_count = 6;

  return (neighbor_count);
}

void grid_changed_box_l() {
  for (int i = 0; i < 3; i++) {
    local_geo.local_box_l[i] = box_geo.length()[i] / (double)node_grid[i];
    local_geo.my_left_[i] = node_pos[i] * local_geo.length()[i];
    local_geo.my_right_[i] = (node_pos[i] + 1) * local_geo.length()[i];
  }

  calc_minimal_box_dimensions();
}

void grid_changed_n_nodes() {
  mpi_reshape_communicator({{node_grid[0], node_grid[1], node_grid[2]}},
                           {{1, 1, 1}});

  MPI_Cart_coords(comm_cart, this_node, 3, node_pos.data());

  calc_node_neighbors(this_node);

  grid_changed_box_l();
}

void calc_minimal_box_dimensions() {
  min_box_l = *boost::min_element(box_geo.length());
  min_local_box_l = *boost::min_element(local_geo.length());
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

void map_node_array(int node, int pos[3]) {
  MPI_Cart_coords(comm_cart, node, 3, pos);
}

int map_array_node(Utils::Span<const int> pos) {
  int rank;
  MPI_Cart_rank(comm_cart, pos.data(), &rank);
  return rank;
}
