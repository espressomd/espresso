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
 *
 * Lattice class definition
 *
 */

#include <boost/range/numeric.hpp>

#include "grid.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include <utils/constants.hpp>
#include <utils/index.hpp>

Lattice::Lattice(double agrid, double offset, int halo_size,
                 Utils::Vector3d const &local_box,
                 Utils::Vector3d const &my_right,
                 Utils::Vector3d const &box_length,
                 Utils::Vector3i const &node_pos,
                 Utils::Vector3i const &node_grid)
    : agrid(agrid), halo_size(halo_size), offset(offset), node_grid(node_grid),
      local_box(local_box), my_right(my_right) {
  /* determine the number of local lattice nodes */
  auto const epsilon = std::numeric_limits<double>::epsilon();
  for (int d = 0; d < 3; d++) {
    grid[d] = static_cast<int>(round(local_box[d] / agrid));
    global_grid[d] = node_grid[d] * grid[d];
    local_index_offset[d] = node_pos[d] * grid[d];
  }

  // sanity checks
  for (int dir = 0; dir < 3; dir++) {
    // check if local_box_l is compatible with lattice spacing
    auto diff = fabs(local_box[dir] - grid[dir] * agrid);
    if (diff > epsilon * box_length[dir]) {
      throw std::runtime_error(
          "Lattice spacing agrid[" + std::to_string(dir) +
          "]=" + std::to_string(agrid) + " is incompatible with local_box_l[" +
          std::to_string(dir) + "]=" + std::to_string(local_box[dir]) +
          " ( box_l[" + std::to_string(dir) +
          "]=" + std::to_string(box_length[dir]) +
          " ). Mismatch: " + std::to_string(diff));
    }
  }

  /* determine the number of total nodes including halo */
  halo_grid = grid + Utils::Vector3i::broadcast(2 * halo_size);
  halo_grid_volume = boost::accumulate(halo_grid, 1, std::multiplies<int>());
  halo_offset =
      Utils::get_linear_index(halo_size, halo_size, halo_size, halo_grid);
}

bool Lattice::is_local(Utils::Vector3i const &index) const noexcept {
  auto const x = index * agrid;
  return x >= my_right - local_box and x < my_right;
}

void Lattice::map_position_to_lattice(const Utils::Vector3d &pos,
                                      Utils::Vector<std::size_t, 8> &node_index,
                                      Utils::Vector6d &delta) const {
  Utils::Vector3i ind{};
  auto const epsilon = std::numeric_limits<double>::epsilon();

  /* determine the elementary lattice cell containing the particle
     and the relative position of the particle in this cell */
  for (int dir = 0; dir < 3; dir++) {
    auto const lpos = pos[dir] - (my_right[dir] - local_box[dir]);
    auto const rel = lpos / agrid + offset;
    ind[dir] = static_cast<int>(floor(rel));

    /* surrounding elementary cell is not completely inside this box,
       adjust if this is due to round off errors */
    if (ind[dir] < 0) {
      if (fabs(rel) < epsilon) {
        ind[dir] = 0; // TODO
      } else {
        throw std::runtime_error("position not inside a local plaquette");
      }
    } else if (ind[dir] > grid[dir]) {
      if (lpos - local_box[dir] < epsilon * local_box[dir])
        ind[dir] = grid[dir];
      else
        throw std::runtime_error("position not inside a local plaquette");
    }

    delta[3 + dir] = rel - ind[dir]; // delta_x/a
    delta[dir] = 1.0 - delta[3 + dir];
  }
  node_index[0] = Utils::get_linear_index(ind, halo_grid);
  node_index[1] = node_index[0] + 1;
  node_index[2] = node_index[0] + halo_grid[0];
  node_index[3] = node_index[0] + halo_grid[0] + 1;
  node_index[4] = node_index[0] + halo_grid[0] * halo_grid[1];
  node_index[5] = node_index[4] + 1;
  node_index[6] = node_index[4] + halo_grid[0];
  node_index[7] = node_index[4] + halo_grid[0] + 1;
}

Utils::Vector3i Lattice::local_index(Utils::Vector3i const &global_index) const
    noexcept {
  return {global_index[0] - local_index_offset[0] + halo_size,
          global_index[1] - local_index_offset[1] + halo_size,
          global_index[2] - local_index_offset[2] + halo_size};
}
