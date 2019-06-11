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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/************************************************
 * defines
 ************************************************/

#define MAX_INTERACTION_RANGE 1e100

/**********************************************
 * variables
 **********************************************/

Utils::Vector3i node_grid{};
Utils::Vector3i node_pos = {-1, -1, -1};
Utils::Vector<int, 6> node_neighbors{};
Utils::Vector<int, 6> boundary{};
int periodic = 7;

Utils::Vector3d box_l = {1, 1, 1};
Utils::Vector3d half_box_l = {0.5, 0.5, 0.5};
Utils::Vector3d box_l_i = {1, 1, 1};
double min_box_l;
Utils::Vector3d local_box_l{1, 1, 1};
double min_local_box_l;
Utils::Vector3d my_left{};
Utils::Vector3d my_right{1, 1, 1};

/************************************************************/

void init_node_grid() {
  grid_changed_n_nodes();
  cells_on_geometry_change(CELL_FLAG_GRIDCHANGED);
}

int map_position_node_array(const Utils::Vector3d &pos) {
  auto const f_pos = folded_position(pos);

  Utils::Vector3i im;
  for (int i = 0; i < 3; i++) {
    im[i] = std::floor(f_pos[i] / local_box_l[i]);
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
      boundary[2 * dir] = 1;
    } else {
      boundary[2 * dir] = 0;
    }
    /* right boundary ? */
    if (node_pos[dir] == node_grid[dir] - 1) {
      boundary[2 * dir + 1] = -1;
    } else {
      boundary[2 * dir + 1] = 0;
    }
  }

  neighbor_count = 6;
  GRID_TRACE(printf("%d: node_grid %d %d %d, pos %d %d %d, node_neighbors ",
                    this_node, node_grid[0], node_grid[1], node_grid[2],
                    node_pos[0], node_pos[1], node_pos[2]));

  return (neighbor_count);
}

void grid_changed_box_l() {
  int i;

  GRID_TRACE(fprintf(stderr, "%d: grid_changed_box_l:\n", this_node));
  GRID_TRACE(fprintf(stderr, "%d: node_pos %d %d %d\n", this_node, node_pos[0],
                     node_pos[1], node_pos[2]));
  GRID_TRACE(fprintf(stderr, "%d: node_grid %d %d %d\n", this_node,
                     node_grid[0], node_grid[1], node_grid[2]));
  for (i = 0; i < 3; i++) {
    local_box_l[i] = box_l[i] / (double)node_grid[i];
    my_left[i] = node_pos[i] * local_box_l[i];
    my_right[i] = (node_pos[i] + 1) * local_box_l[i];
    box_l_i[i] = 1 / box_l[i];
    half_box_l[i] = 0.5 * box_l[i];
  }

  calc_minimal_box_dimensions();

#ifdef GRID_DEBUG
  fprintf(stderr, "%d: local_box_l = (%.3f, %.3f, %.3f)\n", this_node,
          local_box_l[0], local_box_l[1], local_box_l[2]);
  fprintf(stderr,
          "%d: coordinates: x in [%.3f, %.3f], y in [%.3f, %.3f], z in "
          "[%.3f, %.3f]\n",
          this_node, my_left[0], my_right[0], my_left[1], my_right[1],
          my_left[2], my_right[2]);
#endif
}

void grid_changed_n_nodes() {
  GRID_TRACE(fprintf(stderr, "%d: grid_changed_n_nodes:\n", this_node));

  mpi_reshape_communicator({{node_grid[0], node_grid[1], node_grid[2]}},
                           {{1, 1, 1}});

  MPI_Cart_coords(comm_cart, this_node, 3, node_pos.data());

  calc_node_neighbors(this_node);

#ifdef GRID_DEBUG
  fprintf(stderr, "%d: node_pos=(%d,%d,%d)\n", this_node, node_pos[0],
          node_pos[1], node_pos[2]);
  fprintf(stderr, "%d: node_neighbors=(%d,%d,%d,%d,%d,%d)\n", this_node,
          node_neighbors[0], node_neighbors[1], node_neighbors[2],
          node_neighbors[3], node_neighbors[4], node_neighbors[5]);
  fprintf(stderr, "%d: boundary=(%d,%d,%d,%d,%d,%d)\n", this_node, boundary[0],
          boundary[1], boundary[2], boundary[3], boundary[4], boundary[5]);
#endif

  grid_changed_box_l();
}

void calc_minimal_box_dimensions() {
  int i;
  min_box_l = 2 * MAX_INTERACTION_RANGE;
  min_local_box_l = MAX_INTERACTION_RANGE;
  for (i = 0; i < 3; i++) {
    min_box_l = std::min(min_box_l, box_l[i]);
    min_local_box_l = std::min(min_local_box_l, local_box_l[i]);
  }
}

void rescale_boxl(int dir, double d_new) {
  double scale = (dir - 3) ? d_new / box_l[dir] : d_new / box_l[0];
  if (scale < 1.) {
    mpi_rescale_particles(dir, scale);
    if (dir < 3)
      box_l[dir] = d_new;
    else
      box_l[0] = box_l[1] = box_l[2] = d_new;
    mpi_bcast_parameter(FIELD_BOXL);
  } else if (scale > 1.) {
    if (dir < 3)
      box_l[dir] = d_new;
    else
      box_l[0] = box_l[1] = box_l[2] = d_new;
    mpi_bcast_parameter(FIELD_BOXL);
    mpi_rescale_particles(dir, scale);
  }
}
