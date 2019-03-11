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
 * Contains the lattice layout and pointers to the data fields.
 * For parallelization purposes, it is assumed that a halo region
 * surrounds the local lattice sites.
 */

#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include "utils/Vector.hpp"

/** Switch determining the type of lattice dynamics. A value of zero
 *  means that there is no lattice dynamics. Different types can be
 *  combined by or'ing the respective flags.
 *  So far, only \ref LATTICE_OFF and \ref LATTICE_LB exist.
 */
extern int lattice_switch;

#define LATTICE_LB 1     /** Lattice Boltzmann */
#define LATTICE_LB_GPU 2 /** Lattice Boltzmann */
#define LATTICE_OFF 0    /** Lattice off */

class Lattice {
public:
  using index_t = int;

  Vector3i grid; /**< number of local lattice sites in each direction
                  *   (excluding halo) */
  Vector3i global_grid;
  Vector3d agrid; /**< lattice constant */

  Vector3i halo_grid; /**< number of lattice sites in each direction
                       *   (including halo) */
  int halo_size;      /**< halo size in all directions */

  Vector3d offset; /**< global offset */
  Vector3d local_offset;
  Vector3i local_index_offset;

  index_t halo_grid_volume; /**< total number (volume) of lattice sites
                             *   (including halo) */
  index_t halo_offset; /**< offset for number of halo sites stored in front of
                        *   the local lattice sites */

  /** Initialize lattice.
   *
   * This function initializes the variables describing the lattice
   * layout. Important: The lattice data is <em>not</em> allocated here!
   *
   * \param agrid      lattice spacing
   * \param offset     lattice offset
   * \param halo_size  halo size
   * \param dim        lattice dimensions
   */
  int init(double *agrid, double *offset, int halo_size, size_t dim);

  /** Map a spatial position to the surrounding lattice sites.
   *
   * This function takes a global spatial position and determines the
   * surrounding elementary cell of the lattice for this position.
   * The distance fraction in each direction is also calculated.
   * <br><em>Remarks:</em>
   * <ul>
   * <li>The spatial position has to be in the local domain.</li>
   * <li>The lattice sites of the elementary cell are returned as local
   * indices</li>
   * </ul>
   * \param pos        spatial position (Input)
   * \param node_index local indices of the surrounding lattice sites (Output)
   * \param delta      distance fraction of pos from the surrounding
   *                   elementary cell, 6 directions (Output)
   */
  void map_position_to_lattice(const Vector3d &pos,
                               Vector<std::size_t, 8> &node_index,
                               Vector6d &delta) const;
  /********************** Inline Functions **********************/

  /** Map a global lattice site to the node grid.
   *
   *  This function determines the processor responsible for
   *  the specified lattice site. The coordinates of the site are
   *  taken as global coordinates.
   *
   * \param  ind     global coordinates of the lattice site
   * \return         index of the node for the lattice site
   */
  int map_lattice_to_node(Vector3i &ind) const;

  /********************** static Functions **********************/

  /** Map a spatial position to the surrounding lattice sites.
   *
   * This function takes a global spatial position and determines the
   * surrounding elementary cell of the lattice for this position.
   * The distance fraction in each direction is also calculated.
   * <br><em>Remarks:</em>
   * <ul>
   * <li>The spatial position is given in global coordinates.</li>
   * <li>The lattice sites of the elementary cell are returned as local
   * indices</li>
   * </ul>
   * \param pos        spatial position (Input)
   * \param ind        global index of the lower left lattice site (Output)
   * \param delta      distance fraction of pos from the surrounding
   *                   elementary cell, 6 directions (Output)
   * \param tmp_agrid  lattice mesh distance
   */
  static void map_position_to_lattice_global(const Vector3d &pos, Vector3i &ind,
                                             Vector6d &delta, double tmp_agrid);
};

#endif /* LATTICE_HPP */
