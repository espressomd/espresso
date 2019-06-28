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

#include <bitset>
#include <boost/range/numeric.hpp>

#include "grid_based_algorithms/lattice.hpp"

#include "debug.hpp"
#include "grid.hpp"

#include <utils/index.hpp>
using Utils::get_linear_index;
#include <utils/constants.hpp>

bool Lattice::is_local(Utils::Vector3i const &index) const noexcept {
  std::bitset<3> greater_equal_left_bound;
  std::bitset<3> smaller_right_bound;
  double x;
  for (int i = 0; i < index.size(); ++i) {
    x = offset + index[i] * agrid;
    if (x >= my_left[i] + offset) {
      greater_equal_left_bound.set(i);
    }
    if (x < my_right[i] + offset) {
      smaller_right_bound.set(i);
    }
  }
  auto const res =
      (greater_equal_left_bound.all() and smaller_right_bound.all());
  return res;
}

Lattice::Lattice(double agrid, double offset, int halo_size,
                 Utils::Vector3d const &local_box,
                 Utils::Vector3d const &myright,
                 Utils::Vector3d const &box_length,
                 Utils::Vector3i const &node_grid)
    : agrid(agrid), halo_size(halo_size), offset(offset), node_grid(node_grid) {
  /* determine the number of local lattice nodes */
  auto const epsilon = std::numeric_limits<double>::epsilon();
  for (int d = 0; d < 3; d++) {
    global_grid[d] = static_cast<int>(std::round(box_length[d] / agrid));
    local_index_offset[d] =
        static_cast<int>(ceil((my_left[d] - offset) / agrid));
    local_offset[d] = offset + local_index_offset[d] * agrid;
    grid[d] = static_cast<int>(
        ceil((myright[d] - local_offset[d] - epsilon) / agrid));
  }

  // sanity checks
  for (int dir = 0; dir < 3; dir++) {
    // check if local_box_l is compatible with lattice spacing
    if (fabs(local_box[dir] - grid[dir] * agrid) > epsilon * box_length[dir]) {
      throw std::runtime_error(
          "Lattice spacing agrid[" + std::to_string(dir) +
          "]=" + std::to_string(agrid) + " is incompatible with local_box_l[" +
          std::to_string(dir) + "]=" + std::to_string(local_box[dir]) +
          " ( box_l[" + std::to_string(dir) +
          "]=" + std::to_string(box_length[dir]) + " )");
    }
  }

  /* determine the number of total nodes including halo */
  halo_grid = grid + Utils::Vector3i::broadcast(2 * halo_size);

  halo_grid_volume = boost::accumulate(halo_grid, 1, std::multiplies<int>());
  this->halo_offset =
      get_linear_index(halo_size, halo_size, halo_size, this->halo_grid);
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
    } else if (ind[dir] > this->grid[dir]) {
      if (lpos - local_box[dir] < epsilon * local_box[dir])
        ind[dir] = this->grid[dir];
      else
        throw std::runtime_error("position not inside a local plaquette");
    }

    delta[3 + dir] = rel - ind[dir]; // delta_x/a
    delta[dir] = 1.0 - delta[3 + dir];
  }
  node_index[0] = get_linear_index(ind, halo_grid);
  node_index[1] = node_index[0] + 1;
  node_index[2] = node_index[0] + halo_grid[0];
  node_index[3] = node_index[0] + halo_grid[0] + 1;
  node_index[4] = node_index[0] + halo_grid[0] * halo_grid[1];
  node_index[5] = node_index[4] + 1;
  node_index[6] = node_index[4] + halo_grid[0];
  node_index[7] = node_index[4] + halo_grid[0] + 1;
}

int Lattice::map_lattice_to_node(Utils::Vector3i &ind) const noexcept {
  /* determine coordinates in node_grid */
  Utils::Vector3i index_grid;
  index_grid[0] =
      static_cast<int>(floor(ind[0] * agrid / box_l[0] * node_grid[0]));
  index_grid[1] =
      static_cast<int>(floor(ind[1] * agrid / box_l[1] * node_grid[1]));
  index_grid[2] =
      static_cast<int>(floor(ind[2] * agrid / box_l[2] * node_grid[2]));

  /* change from global to local lattice coordinates */
  ind[0] = ind[0] - index_grid[0] * grid[0] + halo_size;
  ind[1] = ind[1] - index_grid[1] * grid[1] + halo_size;
  ind[2] = ind[2] - index_grid[2] * grid[2] + halo_size;

  /* return linear index into node array */
  return map_array_node({grid.data(), 3});
}
