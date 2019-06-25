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
#ifndef _GRID_H
#define _GRID_H
/** \file
 *  Domain decomposition for parallel computing.
 *
 *  The primary simulation box is divided into orthogonal rectangular
 *  subboxes which are assigned to the different nodes (or processes
 *  or threads if you want). This grid is described in \ref
 *  node_grid. Each node has a number \ref this_node and a position
 *  \ref node_pos in that grid. Each node has also 6 nearest neighbors
 *  \ref node_neighbors which are necessary for the communication
 *  between the nodes (see also \ref ghosts.cpp and \ref p3m.cpp for more
 *  details about the communication.
 *
 *  For the 6 directions \anchor directions we have the following convention:
 *
 *  \image html directions.gif "Convention for the order of the directions"
 *
 *  The Figure illustrates the direction convention used for arrays
 *  with 6 (e.g. \ref node_neighbors, \ref #boundary) and 3 entries
 *  (e.g \ref node_grid, box length , \ref my_left,...).
 *
 *
 *  Implementation in grid.cpp.
 */

#include "BoxGeometry.hpp"
#include "algorithm/periodic_fold.hpp"

#include "LocalBox.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <cassert>
#include <limits>

extern BoxGeometry box_geo;
extern LocalBox<double> local_geo;

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** The number of nodes in each spatial dimension. */
extern Utils::Vector3i node_grid;
/** position of node in node grid */
extern Utils::Vector3i node_pos;
/** the six nearest neighbors of a node in the node grid. */
extern Utils::Vector<int, 6> node_neighbors;
/** where to fold particles that leave local box in direction i. */
extern Utils::Vector<int, 6> boundary;

/** Smallest simulation box dimension. Only periodic directions
    are taken into account! */
extern double min_box_l;
/** Dimensions of the box a single node is responsible for. */
extern Utils::Vector3d local_box_l;
/** Smallest local simulation box dimension (\ref local_box_l).
    Only the periodic directions are taken into account! */
extern double min_local_box_l;
/** Right (top, back) corner of this nodes local box. */
extern Utils::Vector3d my_right_global_abcd;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Make sure that the node grid is set, eventually
    determine one automatically. */
void init_node_grid();

/** node mapping: array -> node.
 *
 * \param node   rank of the node you want to know the position for.
 * \param pos    position of the node in node grid.
 */
void map_node_array(int node, int pos[3]);

/** node mapping: node -> array.
 *
 * \return      rank of the node at position pos.
 * \param pos   position of the node in node grid.
 */
int map_array_node(Utils::Span<const int> pos);

/** map a spatial position to the node grid */
int map_position_node_array(const Utils::Vector3d &pos);

/** fill neighbor lists of node.
 *
 * Calculates the numbers of the nearest neighbors for a node and
 * stores them in \ref node_neighbors.
 *
 * \return     the number of neighbors
 * \param node number of the node.  */
int calc_node_neighbors(int node);

/** called from \ref mpi_bcast_parameter . */
void grid_changed_n_nodes();

/** called from \ref mpi_bcast_parameter . */
void grid_changed_box_l();

/** Calculates the smallest box and local box dimensions for periodic
 * directions.  This is needed to check if the interaction ranges are
 * compatible with the box dimensions and the node grid.
 * Remark: In the aperiodic case min_box_l is set to
 * 2 * \ref MAX_INTERACTION_RANGE . */
void calc_minimal_box_dimensions();

/** rescales the box in dimension 'dir' to the new value 'd_new', and rescales
 * the particles accordingly */
void rescale_boxl(int dir, double d_new);

template <typename T> T get_mi_coord(T a, T b, T box_length, bool periodic) {
  auto const dx = a - b;

  if (periodic && (std::fabs(dx) > (0.5 * box_length))) {
    return dx - std::round(dx * (1. / box_length)) * box_length;
  }

  return dx;
}

template <typename T>
Utils::Vector<T, 3> get_mi_vector(const Utils::Vector<T, 3> &a,
                                  const Utils::Vector<T, 3> &b,
                                  const BoxGeometry &box) {
  return {get_mi_coord(a[0], b[0], box.length()[0], box.periodic(0)),
          get_mi_coord(a[1], b[1], box.length()[1], box.periodic(1)),
          get_mi_coord(a[2], b[2], box.length()[2], box.periodic(2))};
}

/** fold a coordinate to primary simulation box.
    \param pos         the position...
    \param image_box   and the box
    \param length the box length.
   coordinate.

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
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

/** fold particle coordinates to primary simulation box.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
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

inline Utils::Vector3d image_shift(const Utils::Vector3i &image_box,
                                   const Utils::Vector3d &box) {
  return {image_box[0] * box[0], image_box[1] * box[1], image_box[2] * box[2]};
}

inline Utils::Vector3d unfolded_position(const Utils::Vector3d &pos,
                                         const Utils::Vector3i &image_box,
                                         const Utils::Vector3d &box) {
  return pos + image_shift(image_box, box);
}

/*@}*/
#endif
