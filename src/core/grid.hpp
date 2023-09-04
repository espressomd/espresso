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

#pragma once

/** @file
 *  Domain decomposition for parallel computing.
 *
 *  The primary simulation box is divided into orthogonal rectangular
 *  subboxes which are assigned to the different MPI nodes using a
 *  Cartesian topolgy.
 *
 *  Implementation in grid.cpp.
 */

#include "BoxGeometry.hpp"

#include "LocalBox.hpp"

#include <utils/Vector.hpp>

extern BoxGeometry box_geo;
extern LocalBox local_geo;

void grid_changed_node_grid(bool update_box_geo = true);

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
 * @param box         Geometry of the simulation box
 * @param node_index  Node index in the Cartesian topology
 * @param node_grid   Dimensions of the Cartesian topology
 * @return Geometry for the node
 */
LocalBox regular_decomposition(BoxGeometry const &box,
                               Utils::Vector3i const &node_index,
                               Utils::Vector3i const &node_grid);

void set_node_grid(Utils::Vector3i const &value);
void set_box_length(Utils::Vector3d const &value);
