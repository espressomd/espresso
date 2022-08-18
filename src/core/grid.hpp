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
#ifndef CORE_GRID_HPP
#define CORE_GRID_HPP
/** @file
 *  Domain decomposition for parallel computing.
 *
 *  The primary simulation box is divided into orthogonal rectangular
 *  subboxes which are assigned to the different nodes (or processes
 *  or threads if you want). This grid is described in @ref node_grid.
 *
 *  Implementation in grid.cpp.
 */

#include "BoxGeometry.hpp"

#include "LocalBox.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

extern BoxGeometry box_geo;
extern LocalBox<double> local_geo;

/** The number of nodes in each spatial dimension. */
extern Utils::Vector3i node_grid;

/** Make sure that the node grid is set, eventually
 *  determine one automatically.
 */
void init_node_grid();

/** @brief Fill neighbor lists of node.
 *
 * Calculates the numbers of the nearest neighbors for a node.
 *
 * @return Ranks of neighbors
 */
Utils::Vector<int, 6> calc_node_neighbors(const boost::mpi::communicator &comm);

/**
 * @brief Calculate the position of node in topology.
 *
 * @param comm Cartesian communicator
 * @return Index of node in grid.
 */
Utils::Vector3i calc_node_pos(const boost::mpi::communicator &comm);

void grid_changed_n_nodes();

void grid_changed_box_l(const BoxGeometry &box);

/** @brief Calculate image box shift vector.
 *  @param image_box  image box offset
 *  @param box        box parameters (side lengths)
 *  @return Image box coordinates.
 */
inline Utils::Vector3d image_shift(const Utils::Vector3i &image_box,
                                   const Utils::Vector3d &box) {
  return hadamard_product(image_box, box);
}

/** @brief Unfold particle coordinates to image box.
 *  @param pos        coordinate to unfold
 *  @param image_box  image box offset
 *  @param box        box parameters (side lengths, periodicity)
 *  @return Unfolded coordinates.
 */
inline Utils::Vector3d unfolded_position(const Utils::Vector3d &pos,
                                         const Utils::Vector3i &image_box,
                                         const Utils::Vector3d &box) {
  return pos + image_shift(image_box, box);
}

/**
 * @brief Composition of the simulation box into equal parts for each node.
 *
 * @param box Geometry of the simulation box
 * @param node_pos Position of node in the node grid
 * @param node_grid Nodes in each direction
 * @return Geometry for the node
 */
LocalBox<double> regular_decomposition(const BoxGeometry &box,
                                       Utils::Vector3i const &node_pos,
                                       Utils::Vector3i const &node_grid);

void set_node_grid(Utils::Vector3i const &value);
void set_box_length(Utils::Vector3d const &value);

#endif
