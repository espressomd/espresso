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
 *  node_grid.
 *
 *  Implementation in grid.cpp.
 */

#include "BoxGeometry.hpp"
#include "algorithm/periodic_fold.hpp"

#include "LocalBox.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/math/sgn.hpp>
#ifdef LEES_EDWARDS
#include "lees_edwards.hpp"
#endif

#include <boost/mpi/communicator.hpp>
#include <cassert>
#include <limits>

extern BoxGeometry box_geo;
extern LocalBox<double> local_geo;

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** The number of nodes in each spatial dimension. */
extern Utils::Vector3i node_grid;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Make sure that the node grid is set, eventually
    determine one automatically. */
void init_node_grid();

/** map a spatial position to the node grid */
int map_position_node_array(const Utils::Vector3d &pos);

/** fill neighbor lists of node.
 *
 * Calculates the numbers of the nearest neighbors for a node.
 *
 * \return Ranks of neighbors
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

/** get the minimal distance vector of two vectors in the current bc.
 *  @param a the vector to subtract from
 *  @param b the vector to subtract
 *  @param res where to store the result
 */

template <typename T, typename U, typename V>
inline void get_mi_vector(T &res, U const &a, V const &b,
                          const BoxGeometry &box) {
#ifdef LEES_EDWARDS
  const double &offset = box.lees_edwards_state.pos_offset;
  const unsigned int shear_plane_normal =
      LeesEdwards::get_shear_plane_normal_coord(box.lees_edwards_protocol);
  const unsigned int shear_dir =
      LeesEdwards::get_shear_dir_coord(box.lees_edwards_protocol);

  const double dist = a[shear_plane_normal] - b[shear_plane_normal];
#endif
  for (int i = 0; i < 3; i++) {
#ifdef LEES_EDWARDS
    double shift;
    if (i == shear_dir &&
        std::fabs(dist) > .5 * box.length()[shear_plane_normal]) {
      shift =
          Utils::sgn(dist) * (offset - round(offset / box.length()[shear_dir]) *
                                           box.length()[shear_dir]);
    } else {
      shift = 0.0;
    }
#else
    constexpr const double shift = 0.;
#endif
    res[i] = get_mi_coord(a[i] - shift, b[i], box.length()[i], box.periodic(i));
  }
}

template <typename T, typename U>
Utils::Vector3d get_mi_vector(T const &a, U const &b, const BoxGeometry &box) {
  Utils::Vector3d res;
  get_mi_vector(res, a, b, box);

  return res;
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
  return hadamard_product(image_box, box);
}

inline Utils::Vector3d unfolded_position(const Utils::Vector3d &pos,
                                         const Utils::Vector3i &image_box,
                                         const Utils::Vector3d &box) {
  return pos + image_shift(image_box, box);
}

/** Calculate the velocity difference including the Lees Edwards velocity*/
inline Utils::Vector3d vel_diff(Utils::Vector3d const &x,
                                Utils::Vector3d const &y,
                                Utils::Vector3d const &u,
                                Utils::Vector3d const &v,
                                const BoxGeometry &box) {

  auto ret = u - v;

#ifdef LEES_EDWARDS
  auto const shear_velocity = box.lees_edwards_state.shear_velocity;
  auto const shear_plane_normal =
      LeesEdwards::get_shear_plane_normal_coord(box.lees_edwards_protocol);
  auto const shear_dir =
      LeesEdwards::get_shear_dir_coord(box.lees_edwards_protocol);
  auto const dy = std::abs(x[shear_plane_normal] - y[shear_plane_normal]);
  if (dy > 0.5 * box.length()[shear_plane_normal]) {
    ret[shear_dir] += Utils::sgn(dy) * shear_velocity;
  }
#endif

  return ret;
}

class PositionFolder {
public:
  template <typename Particle> void operator()(Particle &p) const {
    fold_position(p.r.p, p.l.i);
  }
};

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
/*@}*/
#endif
