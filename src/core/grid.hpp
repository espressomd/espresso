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
#ifndef _GRID_H
#define _GRID_H
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
#include "algorithm/periodic_fold.hpp"

#include "LocalBox.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>
#include <cassert>
#include <limits>

extern BoxGeometry box_geo;
extern LocalBox<double> local_geo;

/** The number of nodes in each spatial dimension. */
extern Utils::Vector3i node_grid;

/** Make sure that the node grid is set, eventually
 *  determine one automatically.
 */
void init_node_grid();

/** @brief Map a spatial position to the node grid */
int map_position_node_array(const Utils::Vector3d &pos);

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

/** called from \ref mpi_bcast_parameter . */
void grid_changed_n_nodes();

/** called from \ref mpi_bcast_parameter . */
void grid_changed_box_l(const BoxGeometry &box);

/** @brief Rescale box in dimension @p dir to the new value @p d_new and
 *  rescale the particles accordingly.
 */
void rescale_boxl(int dir, double d_new);

/**
 * @brief Get the minimum-image distance between two coordinates.
 * @param a           Coordinate of the terminal point.
 * @param b           Coordinate of the initial point.
 * @param box_length  Box length.
 * @param periodic    Box periodicity.
 * @return Shortest distance from @p b to @p a across periodic images,
 *         i.e. <tt>a - b</tt>. Can be negative.
 */
template <typename T> T get_mi_coord(T a, T b, T box_length, bool periodic) {
  auto const dx = a - b;

  if (periodic && (std::fabs(dx) > (0.5 * box_length))) {
    return dx - std::round(dx * (1. / box_length)) * box_length;
  }

  return dx;
}

/**
 * @brief Get the minimum-image vector between two coordinates.
 * @param a     Coordinate of the terminal point.
 * @param b     Coordinate of the initial point.
 * @param box   Box parameters (side lengths, periodicity).
 * @return Vector from @p b to @p a that minimizes the distance across
 *         periodic images, i.e. <tt>a - b</tt>.
 */
template <typename T>
Utils::Vector<T, 3> get_mi_vector(const Utils::Vector<T, 3> &a,
                                  const Utils::Vector<T, 3> &b,
                                  const BoxGeometry &box) {
  return {get_mi_coord(a[0], b[0], box.length()[0], box.periodic(0)),
          get_mi_coord(a[1], b[1], box.length()[1], box.periodic(1)),
          get_mi_coord(a[2], b[2], box.length()[2], box.periodic(2))};
}

/** @brief Fold a coordinate to primary simulation box.
 *  @param pos        coordinate to fold
 *  @param image_box  image box offset
 *  @param length     box length
 */
inline std::pair<double, int> fold_coordinate(double pos, int image_box,
                                              double const &length) {
  std::tie(pos, image_box) = Algorithm::periodic_fold(pos, image_box, length);

  if ((image_box == std::numeric_limits<int>::min()) ||
      (image_box == std::numeric_limits<int>::max())) {
    throw std::runtime_error(
        "Overflow in the image box count while folding a particle coordinate "
        "into the primary simulation box. Maybe a particle experienced a "
        "huge force.");
  }

  return {pos, image_box};
}

/** @brief Fold particle coordinates to primary simulation box.
 *  @param[in,out] pos        coordinate to fold
 *  @param[in,out] image_box  image box offset
 *  @param[in] box            box parameters (side lengths, periodicity)
 */
inline void fold_position(Utils::Vector3d &pos, Utils::Vector3i &image_box,
                          const BoxGeometry &box) {
  for (int i = 0; i < 3; i++) {
    if (box.periodic(i)) {
      std::tie(pos[i], image_box[i]) =
          fold_coordinate(pos[i], image_box[i], box.length()[i]);
    }
  }
}

/** @brief Fold particle coordinates to primary simulation box.
 *  @param p    coordinate to fold
 *  @param box  box parameters (side lengths, periodicity)
 *  @return Folded coordinates.
 */
inline Utils::Vector3d folded_position(const Utils::Vector3d &p,
                                       const BoxGeometry &box) {
  Utils::Vector3d p_folded;
  for (int i = 0; i < 3; i++) {
    if (box.periodic(i)) {
      p_folded[i] = Algorithm::periodic_fold(p[i], box.length()[i]);
    } else {
      p_folded[i] = p[i];
    }
  }

  return p_folded;
}

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
#endif
