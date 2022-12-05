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
 *
 * Lattice class definition.
 * Contains the lattice layout.
 * For parallelization purposes, it is assumed that a halo region
 * surrounds the local lattice sites.
 */

#ifndef CORE_LB_LATTICE_HPP
#define CORE_LB_LATTICE_HPP

#include <utils/Vector.hpp>

#include <cstddef>

class Lattice {
public:
  using index_t = int;

  /** number of local lattice sites in each direction (excluding halo) */
  Utils::Vector3i grid;
  Utils::Vector3i global_grid;
  double agrid; /**< lattice constant */

  /** number of lattice sites in each direction (including halo) */
  Utils::Vector3i halo_grid;
  index_t halo_size; /**< halo size in all directions */

  double offset; /**< global offset */
  /** global index of the local domain origin */
  Utils::Vector3i local_index_offset;
  /** global domain partition */
  Utils::Vector3i node_grid;
  /** dimensions of the local domain */
  Utils::Vector3d local_box;
  /** global position of the top right corner of the local domain */
  Utils::Vector3d my_right;

  /** total number of lattice sites (including halo) */
  index_t halo_grid_volume;
  /** offset for number of halo sites stored in front of the local
   *  lattice sites
   */
  index_t halo_offset;

  Lattice() = default;
  /** @brief %Lattice constructor.
   *
   *  This function initializes the variables describing the lattice
   *  layout. Important: The lattice data is <em>not</em> allocated here!
   *
   *  @param agrid       lattice spacing
   *  @param offset      lattice offset
   *  @param halo_size   halo size
   *  @param local_box   dimensions of the local box
   *  @param myright     right (top, back) corner of the local box
   *  @param box_length  lengths of the local box
   *  @param node_pos    position of this node in the domain decomposition
   *  @param node_grid   node_grid of domain decomposition
   */
  Lattice(double agrid, double offset, int halo_size,
          const Utils::Vector3d &local_box, const Utils::Vector3d &myright,
          const Utils::Vector3d &box_length, Utils::Vector3i const &node_pos,
          Utils::Vector3i const &node_grid);

  /** Map a spatial position to the surrounding lattice sites.
   *
   * This function takes a global spatial position and determines the
   * surrounding elementary cell of the lattice for this position.
   * The distance fraction in each direction is also calculated.
   *
   * Remarks:
   * - The spatial position has to be in the local domain
   * - The lattice sites of the elementary cell are returned as local indices
   *
   * @param[in]  pos        spatial position
   * @param[out] node_index local indices of the surrounding lattice sites
   * @param[out] delta      distance fraction of %p pos from the surrounding
   *                        elementary cell, 6 directions
   */
  void map_position_to_lattice(Utils::Vector3d const &pos,
                               Utils::Vector<std::size_t, 8> &node_index,
                               Utils::Vector6d &delta) const;

  /**
   * @brief Determine if given global index is node-local.
   * @param index Global lattice index.
   */
  bool is_local(Utils::Vector3i const &index) const noexcept;
  /**
   * @brief Calculate the node-local index.
   * @param global_index Index into global lattice.
   */
  Utils::Vector3i
  local_index(Utils::Vector3i const &global_index) const noexcept;
};

#endif /* CORE_LB_LATTICE_HPP */
