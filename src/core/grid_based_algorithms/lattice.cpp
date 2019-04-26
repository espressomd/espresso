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
 *
 * Lattice class definition
 *
 */

#include "grid_based_algorithms/lattice.hpp"

#include "debug.hpp"
#include "grid.hpp"

#include <utils/index.hpp>
using Utils::get_linear_index;
#include <utils/constants.hpp>

int Lattice::init(double *agrid, double const *offset, int halo_size,
                  const Utils::Vector3d &local_box,
                  const Utils::Vector3d &myright,
                  const Utils::Vector3d &box_length) {
  /* determine the number of local lattice nodes */
  auto const epsilon = std::numeric_limits<double>::epsilon();
  for (int d = 0; d < 3; d++) {
    this->agrid[d] = agrid[d];
    this->global_grid[d] = (int)std::round(box_length[d] / agrid[d]);
    this->offset[d] = offset[d];
    this->local_index_offset[d] =
        (int)ceil((my_left[d] - this->offset[d]) / this->agrid[d]);
    this->local_offset[d] =
        this->offset[d] + this->local_index_offset[d] * this->agrid[d];
    this->grid[d] = (int)ceil((myright[d] - this->local_offset[d] - epsilon) /
                              this->agrid[d]);
  }

  // sanity checks
  for (int dir = 0; dir < 3; dir++) {
    // check if local_box_l is compatible with lattice spacing
    if (fabs(local_box[dir] - this->grid[dir] * agrid[dir]) >
        epsilon * box_length[dir]) {
      runtimeErrorMsg() << "Lattice spacing agrid[" << dir << "]=" << agrid[dir]
                        << " is incompatible with local_box_l[" << dir
                        << "]=" << local_box[dir] << " ( box_l[" << dir
                        << "]=" << box_length[dir] << " )";
      return ES_ERROR;
    }
  }

  LATTICE_TRACE(fprintf(stderr,
                        "%d: box_l (%.3f,%.3f,%.3f) grid (%d,%d,%d) "
                        "node_neighbors (%d,%d,%d,%d,%d,%d)\n",
                        this_node, local_box[0], local_box[1], local_box[2],
                        this->grid[0], this->grid[1], this->grid[2],
                        node_neighbors[0], node_neighbors[1], node_neighbors[2],
                        node_neighbors[3], node_neighbors[4],
                        node_neighbors[5]));

  this->halo_size = halo_size;
  /* determine the number of total nodes including halo */
  this->halo_grid[0] = this->grid[0] + 2 * halo_size;
  this->halo_grid[1] = this->grid[1] + 2 * halo_size;
  this->halo_grid[2] = this->grid[2] + 2 * halo_size;

  this->halo_grid_volume =
      this->halo_grid[0] * this->halo_grid[1] * this->halo_grid[2];
  this->halo_offset =
      get_linear_index(halo_size, halo_size, halo_size, this->halo_grid);

  return ES_OK;
}

void Lattice::map_position_to_lattice(const Utils::Vector3d &pos,
                                      Utils::Vector<std::size_t, 8> &node_index,
                                      Utils::Vector6d &delta,
                                      const Utils::Vector3d &myLeft,
                                      const Utils::Vector3d &local_box) const {
  Utils::Vector3i ind{};
  auto const epsilon = std::numeric_limits<double>::epsilon();

  /* determine the elementary lattice cell containing the particle
     and the relative position of the particle in this cell */
  for (int dir = 0; dir < 3; dir++) {
    auto const lpos = pos[dir] - myLeft[dir];
    auto const rel = lpos / this->agrid[dir] + 0.5; // +1 for halo offset
    ind[dir] = static_cast<int>(floor(rel));

    /* surrounding elementary cell is not completely inside this box,
       adjust if this is due to round off errors */
    if (ind[dir] < 0) {
      if (fabs(rel) < epsilon) {
        ind[dir] = 0; // TODO
      } else {
        fprintf(stderr,
                "%d: map_position_to_lattice: position (%f,%f,%f) not inside a "
                "local plaquette in dir %d ind[dir]=%d rel=%f lpos=%f.\n",
                this_node, pos[0], pos[1], pos[2], dir, ind[dir], rel, lpos);
      }
    } else if (ind[dir] > this->grid[dir]) {
      if (lpos - local_box[dir] < epsilon * local_box[dir])
        ind[dir] = this->grid[dir];
      else
        fprintf(stderr,
                "%d: map_position_to_lattice: position (%f,%f,%f) not inside a "
                "local plaquette in dir %d ind[dir]=%d rel=%f lpos=%f.\n",
                this_node, pos[0], pos[1], pos[2], dir, ind[dir], rel, lpos);
    }

    delta[3 + dir] = rel - ind[dir]; // delta_x/a
    delta[dir] = 1.0 - delta[3 + dir];
  }
  node_index[0] = get_linear_index(ind, this->halo_grid);
  node_index[1] = node_index[0] + 1;
  node_index[2] = node_index[0] + this->halo_grid[0];
  node_index[3] = node_index[0] + this->halo_grid[0] + 1;
  node_index[4] = node_index[0] + this->halo_grid[0] * this->halo_grid[1];
  node_index[5] = node_index[4] + 1;
  node_index[6] = node_index[4] + this->halo_grid[0];
  node_index[7] = node_index[4] + this->halo_grid[0] + 1;
}

int Lattice::map_lattice_to_node(Utils::Vector3i &ind,
                                 const Utils::Vector3i &local_node_grid) const {
  /* determine coordinates in node_grid */
  Utils::Vector3i grid;
  grid[0] =
      (int)floor(ind[0] * this->agrid[0] * box_l_i[0] * local_node_grid[0]);
  grid[1] =
      (int)floor(ind[1] * this->agrid[1] * box_l_i[1] * local_node_grid[1]);
  grid[2] =
      (int)floor(ind[2] * this->agrid[2] * box_l_i[2] * local_node_grid[2]);

  /* change from global to local lattice coordinates */
  ind[0] = ind[0] - grid[0] * this->grid[0] + this->halo_size;
  ind[1] = ind[1] - grid[1] * this->grid[1] + this->halo_size;
  ind[2] = ind[2] - grid[2] * this->grid[2] + this->halo_size;

  /* return linear index into node array */
  return map_array_node({grid.data(), 3});
}
