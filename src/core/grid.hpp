/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file grid.hpp   Domain decomposition for parallel computing.
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
 *  The Figure illustrates the direction convetion used for arrays
 *  with 6 (e.g. \ref node_neighbors, \ref #boundary) and 3 entries
 *  (e.g \ref node_grid, \ref box_l , \ref my_left,...).
 *
 *
 *  For more information on the domain decomposition, see \ref grid.cpp
 * "grid.c".
*/
#include "RuntimeErrorStream.hpp"
#include "communication.hpp"
#include "utils.hpp"
#include "errorhandling.hpp"
#include "lees_edwards.hpp"

#include <climits>

/** Macro that tests for a coordinate being periodic or not. */
#ifdef PARTIAL_PERIODIC
#define PERIODIC(coord) (periodic & (1L << coord))
#else
#define PERIODIC(coord) 1
#endif

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** The number of nodes in each spatial dimension. */
extern int node_grid[3];
/** position of node in node grid */
extern int node_pos[3];
#ifdef LEES_EDWARDS
/** the nearest neighbors of a node in the node grid. */
extern int *node_neighbors;
extern int *node_neighbor_lr;
extern int *node_neighbor_wrap;
/** the number of nearest neighbors of a node in the node grid. */
extern int my_neighbor_count;
#else
/** the six nearest neighbors of a node in the node grid. */
extern int node_neighbors[6];
#endif
/** where to fold particles that leave local box in direction i. */
extern int boundary[6];
/** Flags for all three dimensions wether pbc are applied (default).
    The first three bits give the periodicity */
extern int periodic;

/** Simulation box dimensions. */
extern double box_l[3];
/** Half the box dimensions. Used for get_mi_vector. */
extern double half_box_l[3];
/** 1 / box dimensions. */
extern double box_l_i[3];
/** Smallest simulation box dimension (\ref box_l).
    Remark: with PARTIAL_PERIODIC, only the periodic directions
    are taken into account! */
extern double min_box_l;
/** Dimensions of the box a single node is responsible for. */
extern double local_box_l[3];
/** Smallest local simulation box dimension (\ref local_box_l).
    Remark: with PARTIAL_PERIODIC, only the periodic directions
    are taken into account! */
extern double min_local_box_l;
/** Left (bottom, front) corner of this nodes local box. */
extern double my_left[3];
/** Right (top, back) corner of this nodes local box. */
extern double my_right[3];

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Make sure that the node grid is set, eventually
    determine one automatically. */
void init_node_grid();

/** return wether node grid was set. */
int node_grid_is_set();

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
inline int map_array_node(int pos[3]) {
  int rank;
  MPI_Cart_rank(comm_cart, pos, &rank);
  return rank;
}

/** map a spatial position to the node grid */
int map_position_node_array(double pos[3]);

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
 * Remark: In the apreiodic case min_box_l is set to
 * 2 * \ref MAX_INTERACTION_RANGE . */
void calc_minimal_box_dimensions();

/** calculate most square 2d grid. */
void calc_2d_grid(int n, int grid[3]);

/** calculate 'best' mapping between a 2d and 3d grid.
 *  This we need for the communication from 3d domain decomposition
 *  to 2d row decomposition.
 *  The dimensions of the 2d grid are resorted, if necessary, in a way
 *  that they are multiples of the 3d grid dimensions.
 *  \param g3d      3d grid.
 *  \param g2d      2d grid.
 *  \param mult     factors between 3d and 2d grid dimensions
 *  \return         index of the row direction [0,1,2].
*/
int map_3don2d_grid(int g3d[3], int g2d[3], int mult[3]);

/** rescales the box in dimension 'dir' to the new value 'd_new', and rescales
 * the particles accordingly */
void rescale_boxl(int dir, double d_new);

/** get the minimal distance vector of two vectors in the current bc.
  *  \ref LEES_EDWARDS note: there is no need to add the le_offset here,
  *  any offset should already have been added when the image particle was
 * prepared.
  *  @param a the vector to subtract from
  *  @param b the vector to subtract
  *  @param res where to store the result
*/

template <typename T, typename U, typename V>
inline void get_mi_vector(T &res, U const &a, V const &b) {
  for (int i = 0; i < 3; i++) {
    res[i] = a[i] - b[i];
    if (std::fabs(res[i]) > half_box_l[i] && PERIODIC(i))
      res[i] -= dround(res[i] * box_l_i[i]) * box_l[i];
  }
}

template <typename T, typename U>
Vector3d get_mi_vector(T const &a, U const &b) {
  Vector3d res;
  get_mi_vector(res, a, b);

  return res;
}

/** fold a coordinate to primary simulation box, including velocity (in case of
   LEES_EDWARDS,
    the velocity depends on the image we are in).
    \param pos         the position...
    \param vel         the velocity...
    \param image_box   and the box
    \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z
   coordinate.

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
inline void fold_coordinate(double pos[3], double vel[3], int image_box[3],
                            int dir) {
  if (PERIODIC(dir)) {
    int img_count = (int)floor(pos[dir] * box_l_i[dir]);
    image_box[dir] += img_count;
    pos[dir] = pos[dir] - img_count * box_l[dir];

    if (pos[dir] * box_l_i[dir] < -ROUND_ERROR_PREC ||
        pos[dir] * box_l_i[dir] >= 1 + ROUND_ERROR_PREC) {

      runtimeErrorMsg() << "particle coordinate out of range, pos = "
                        << pos[dir] << ", image box = " << image_box[dir];

      image_box[dir] = 0;
      pos[dir] = 0;
      return;
    }

#ifdef LEES_EDWARDS
    if (dir == 1) {
      /* must image y and v_x at same time as x */
      pos[0] -= (lees_edwards_offset * img_count);
      vel[0] -= (lees_edwards_rate * img_count);

      /* (re)-image x */
      img_count = (int)floor(pos[0] * box_l_i[0]);
      image_box[0] += img_count;
      pos[0] = pos[0] - img_count * box_l[0];
    }

#endif
  }
}

/** fold a coordinate to primary simulation box, not caring about velocities.
    \param pos         the position...
    \param image_box   and the box
    \param dir         the coordinate to fold: dir = 0,1,2 for x, y and z
   coordinate.

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
inline void fold_coordinate(double pos[3], int image_box[3], int dir) {
  double v[3];
  fold_coordinate(pos, v, image_box, dir);
}

/** fold particle coordinates to primary simulation box.
    \param pos the position...
    \param vel the velocity...
    \param image_box and the box

    Pos, vel and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
inline void fold_position(double pos[3], double vel[3], int image_box[3]) {
  for (int i = 0; i < 3; i++)
    fold_coordinate(pos, vel, image_box, i);
}

/** fold particle coordinates to primary simulation box.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O,
    i. e. a previously folded position will be folded correctly.
*/
inline void fold_position(double pos[3], int image_box[3]) {
  for (int i = 0; i < 3; i++)
    fold_coordinate(pos, image_box, i);
}

/** fold particle coordinates to primary simulation box.
 * The particle is not changed.
 */
template <typename Particle> Vector3d folded_position(Particle const &p) {
  Vector3d pos{p.r.p};

  for (int dir = 0; dir < 3; dir++) {
    if (PERIODIC(dir)) {
      const int img_count =
          static_cast<int>(std::floor(pos[dir] * box_l_i[dir]));

      pos[dir] -= img_count * box_l[dir];
    }
  }

  return pos;
}

/** @overload */
template <typename Particle> Vector3d folded_position(const Particle *p) {
  assert(p);
  return folded_position(*p);
}

inline void fold_position(Vector3d &pos, Vector<3, int> &image_box) {
  fold_position(pos.data(), image_box.data());
}

/** unfold coordinates to physical position.
    \param pos the position
    \param pos the velocity
    \param image_box and the box

    Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
    afterwards.
*/
inline void unfold_position(double pos[3], double vel[3], int image_box[3]) {
#ifdef LEES_EDWARDS
  auto const y_img_count = static_cast<int>(floor(pos[1] * box_l_i[1] + image_box[1]));

  pos[0] += image_box[0] * box_l[0] + y_img_count * lees_edwards_offset;
  pos[1] += image_box[1] * box_l[1];
  pos[2] += image_box[2] * box_l[2];

  vel[0] += y_img_count * lees_edwards_rate;

  image_box[0] = image_box[1] = image_box[2] = 0;

#else

  int i;
  for (i = 0; i < 3; i++) {
    pos[i] = pos[i] + image_box[i] * box_l[i];
    image_box[i] = 0;
  }

#endif
}

template<typename Particle>
Vector3d unfolded_position(const Particle * p) {
  Vector3d pos{p->r.p};
#ifdef LEES_EDWARDS
  auto const y_img_count = static_cast<int>(floor(pos[1] * box_l_i[1] + p->l.i[1]));

  pos[0] += p->l.i[0] * box_l[0] + y_img_count * lees_edwards_offset;
  pos[1] += p->l.i[1] * box_l[1];
  pos[2] += p->l.i[2] * box_l[2];
#else
  for (int i = 0; i < 3; i++) {
    pos[i] += p->l.i[i] * box_l[i];
  }

#endif

  return pos;
}

/** unfold coordinates to physical position.
    \param pos the position...
    \param image_box and the box

    Both pos and image_box are I/O, i.e. image_box will be (0,0,0)
    afterwards.
*/
inline void unfold_position(double pos[3], int image_box[3]) {
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
