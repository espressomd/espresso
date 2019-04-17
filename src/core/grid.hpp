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
 *  (e.g \ref node_grid, \ref box_l , \ref my_left,...).
 *
 *
 *  Implementation in grid.cpp.
 */

#include "RuntimeErrorStream.hpp"
#include "algorithm/periodic_fold.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "utils.hpp"
#include "utils/Span.hpp"
#include "utils/Vector.hpp"

#include <limits>

/** Macro that tests for a coordinate being periodic or not. */
#ifdef PARTIAL_PERIODIC
#define PERIODIC(coord) (periodic & (1L << (coord)))
#else
#define PERIODIC(coord) (1)
#endif

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
/** Flags for all three dimensions whether pbc are applied (default).
    The first three bits give the periodicity */
extern int periodic;

/** Simulation box dimensions. */
extern Utils::Vector3d box_l;
/** Half the box dimensions. Used for get_mi_vector. */
extern Utils::Vector3d half_box_l;
/** 1 / box dimensions. */
extern Utils::Vector3d box_l_i;
/** Smallest simulation box dimension (\ref box_l).
    Remark: with PARTIAL_PERIODIC, only the periodic directions
    are taken into account! */
extern double min_box_l;
/** Dimensions of the box a single node is responsible for. */
extern Utils::Vector3d local_box_l;
/** Smallest local simulation box dimension (\ref local_box_l).
    Remark: with PARTIAL_PERIODIC, only the periodic directions
    are taken into account! */
extern double min_local_box_l;
/** Left (bottom, front) corner of this nodes local box. */
extern Utils::Vector3d my_left;
/** Right (top, back) corner of this nodes local box. */
extern Utils::Vector3d my_right;

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
inline void map_node_array(int node, int pos[3]) {
  MPI_Cart_coords(comm_cart, node, 3, pos);
}

/** node mapping: node -> array.
 *
 * \return      rank of the node at position pos.
 * \param pos   position of the node in node grid.
 */
inline int map_array_node(Utils::Span<const int> pos) {
  int rank;
  MPI_Cart_rank(comm_cart, pos.data(), &rank);
  return rank;
}

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
 * see also \ref box_l, \ref local_box_l, \ref min_box_l
 * and \ref min_local_box_l.
 * Remark: In the aperiodic case min_box_l is set to
 * 2 * \ref MAX_INTERACTION_RANGE . */
void calc_minimal_box_dimensions();

/** rescales the box in dimension 'dir' to the new value 'd_new', and rescales
 * the particles accordingly */
void rescale_boxl(int dir, double d_new);

template <typename T> T get_mi_coord(T a, T b, int dir) {
  auto dx = a - b;

  if (PERIODIC(dir) && std::fabs(dx) > half_box_l[dir])
    dx -= std::round(dx * box_l_i[dir]) * box_l[dir];

  return dx;
}

/** get the minimal distance vector of two vectors in the current bc.
 *  @param a the vector to subtract from
 *  @param b the vector to subtract
 *  @param res where to store the result
 */

template <typename T, typename U, typename V>
inline void get_mi_vector(T &res, U const &a, V const &b) {
  for (int i = 0; i < 3; i++) {
    res[i] = get_mi_coord(a[i], b[i], i);
  }
}

template <typename T, typename U>
Utils::Vector3d get_mi_vector(T const &a, U const &b) {
  Utils::Vector3d res;
  get_mi_vector(res, a, b);

  return res;
}

/** fold a coordinate to primary simulation box.
    \param pos         the position...
    \param image_box   and the box
    \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z
   coordinate.

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
template <size_t N, typename T1, typename T2>
void fold_coordinate(Utils::Vector<T1, N> &pos, Utils::Vector<T2, N> &image_box,
                     int dir) {
  if (PERIODIC(dir)) {
    std::tie(pos[dir], image_box[dir]) =
        Algorithm::periodic_fold(pos[dir], image_box[dir], box_l[dir]);

    if ((image_box[dir] == std::numeric_limits<T2>::min()) ||
        (image_box[dir] == std::numeric_limits<T2>::max())) {
      throw std::runtime_error(
          "Overflow in the image box count while folding a particle coordinate "
          "into the primary simulation box. Maybe a particle experienced a "
          "huge force.");
    }
  }
}

/** fold particle coordinates to primary simulation box.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
template <size_t N, typename T1, typename T2>
void fold_position(Utils::Vector<T1, N> &pos, Utils::Vector<T2, N> &image_box) {
  for (int i = 0; i < 3; i++)
    fold_coordinate(pos, image_box, i);
}

inline Utils::Vector3d folded_position(const Utils::Vector3d &p) {
  Utils::Vector3d p_folded;
  for (int i = 0; i < 3; i++) {
    if (PERIODIC(i)) {
      p_folded[i] = Algorithm::periodic_fold(p[i], box_l[i]);
    } else {
      p_folded[i] = p[i];
    }
  }

  return p_folded;
}

/** fold particle coordinates to primary simulation box.
 * The particle is not changed.
 */
inline Utils::Vector3d folded_position(Particle const &p) {
  return folded_position(p.r.p);
}

/** @overload */
inline Utils::Vector3d folded_position(const Particle *p) {
  assert(p);
  return folded_position(p->r.p);
}

/** unfold coordinates to physical position.
    \param pos the position
    \param vel the velocity
    \param image_box and the box

    Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
    afterwards.
*/
template <typename T1, typename T2, typename T3>
void unfold_position(T1 &pos, T2 &vel, T3 &image_box) {
  int i;
  for (i = 0; i < 3; i++) {
    pos[i] = pos[i] + image_box[i] * box_l[i];
    image_box[i] = 0;
  }
}

inline Utils::Vector3d unfolded_position(Particle const *p) {
  Utils::Vector3d pos{p->r.p};
  for (int i = 0; i < 3; i++) {
    pos[i] += p->l.i[i] * box_l[i];
  }

  return pos;
}

inline Utils::Vector3d unfolded_position(Particle const &p) {
  return unfolded_position(&p);
}

/** unfold coordinates to physical position.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
    afterwards.
*/
template <typename T1, typename T2>
void unfold_position(T1 &pos, T2 &image_box) {
  double v[3];
  unfold_position(pos, v, image_box);
}

class PositionUnfolder {
public:
  template <typename Particle> void operator()(Particle &p) const {
    unfold_position(p.r.p, p.l.i);
  }
};

class PositionFolder {
public:
  template <typename Particle> void operator()(Particle &p) const {
    fold_position(p.r.p, p.l.i);
  }
};

/*@}*/
#endif
