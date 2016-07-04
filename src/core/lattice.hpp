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
/** \file lattice.hpp 
 *
 * Lattice class definition
 * Contains the lattice layout and pointers to the data fields.
 * For parallelization purposes, it is assumed that a halo region
 * surrounds the local lattice sites.
 */

#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include "grid.hpp"
#include "particle_data.hpp"

/** Switch determining the type of lattice dynamics. A value of zero
 *  means that there is no lattice dynamics. Different types can be
 *  combined by or'ing the respective flags.
 *  So far, only \ref LATTICE_OFF and \ref LATTICE_LB exist.
 */
extern int lattice_switch;

#define LATTICE_LB   1  /** Lattice Boltzmann */
#define LATTICE_LB_GPU   2  /** Lattice Boltzmann */
#define INTERPOLATION_LINEAR 1
#define index_t long
#define LATTICE_OFF  0  /** Lattice off */

enum { LATTICE_ANISOTROPIC = 1, 
       LATTICE_X_NOTEXT = 2, LATTICE_Y_NOTEXT = 4, LATTICE_Z_NOTEXT = 8 };

class Lattice {
public:
    int grid[3] ;/** number of local lattice sites in each direction (excluding halo) */
    int global_grid[3];
    unsigned int dim;
    double agrid[3];/** lattice constant */

    int halo_grid[3] ;/** number of lattice sites in each direction including halo */
    int halo_size;/** halo size in all directions */

    double offset[3];/** global offset */
    double local_offset[3];
    int local_index_offset[3];

    unsigned int interpolation_type;
    char flags;
    size_t element_size;/** Size of each element in size units (=bytes) */
    size_t lattice_dim;/** Dimension of the field, assuming entries are arrays */

    index_t grid_volume;/** total number (volume) of local lattice sites (excluding halo) */
    index_t halo_grid_volume;/** total number (volume) of lattice sites (including halo) */
    index_t halo_grid_surface;/** number of lattice sites in the halo region */
    index_t halo_offset;/** offset for number of halo sites stored in front of the local lattice sites */

    void *_data;/** pointer to the actual lattice data. This can be a contiguous field of arbitrary data. */

    /** particle representation of this lattice. This is needed to
     *  specify interactions between particles and the lattice.
     *  Actually used are only the identity and the type. */
    Particle part_rep;

    /* Constructor */
    Lattice() {}

    /** Initialize lattice.
     *
     * This function initializes the variables describing the lattice
     * layout. Important: The lattice data is <em>not</em> allocated here!
     *
     * \param lattice pointer to the lattice
     * \param agrid   lattice spacing
     */
    int init(double* agrid, double* offset, int halo_size, size_t dim);

    /** lattice memory allocation.
     * \param lattice pointer to the lattice
     */
    void allocate_memory();
    void allocate_memory(size_t element_size);

    void interpolate(double* pos, double* value);

    void interpolate_linear(double* pos, double* value);

    void interpolate_gradient(double* pos, double* value);

    void interpolate_linear_gradient(double* pos, double* value);

    void set_data_for_global_position_with_periodic_image(double* pos, void* data);

    int global_pos_in_local_box(double pos[3]);

    int global_pos_in_local_halobox(double pos[3]);

    int global_pos_to_lattice_index_checked(double pos[3], int* index);

    /** Map a local lattice site to the global position.
     *
     *  This function determines the processor responsible for
     *  the specified lattice site. The coordinates of the site are
     *  taken as global coordinates andcoordinates of the site are
     *  taken as global coordinates and are returned as local coordinates.
     *
     * \param  lattice pointer to the lattice
     * \param  ind     global coordinates of the lattice site (Input)
     * \param  grid    local coordinates of the lattice site (Output)
     * \return         index of the node for the lattice site
     */
    /* @TODO: Implement! */
    int map_lattice_to_position(int *ind, int *grid);

    void map_linear_index_to_global_pos(index_t ind, double* pos);

    void map_local_index_to_pos(index_t* index, double* pos);

    int map_global_index_to_halo_index(index_t* global_index, index_t* halo_index);

    void map_halo_index_to_pos(index_t* index_in_halogrid, double* pos);

    /** Map a spatial position to the surrounding lattice sites.
     *
     * This function takes a global spatial position and determines the
     * surrounding elementary cell of the lattice for this position.
     * The distance fraction in each direction is also calculated.
     * <br><em>Remarks:</em>
     * <ul>
     * <li>The spatial position has to be in the local domain.</li>
     * <li>The lattice sites of the elementary cell are returned as local indices</li>
     * </ul>
     * \param lattice    pointer to the lattice (Input)
     * \param pos        spatial position (Input)
     * \param node_index local indices of the surrounding lattice sites (Output)
     * \param delta      distance fraction of pos from the surrounding
     *                   elementary cell, 6 directions (Output)
     */
    void map_position_to_lattice(const double pos[3], index_t node_index[8], double delta[6]);

    void get_data_for_halo_index(index_t* ind, void** data);

    void get_data_for_linear_index(index_t ind, void** data);

    void get_data_for_local_index(index_t* ind, void** data);

    void set_data_for_local_halo_grid_index(index_t* ind, void* data);

    void set_data_for_local_grid_index(index_t* ind, void* data);

    int global_pos_to_lattice_halo_index(double* pos, index_t*  ind);

    /********************** Inline Functions **********************/

    /** Map a global lattice site to the node grid.
     *
     *  This function determines the processor responsible for
     *  the specified lattice site. The coordinates of the site are
     *  taken as global coordinates and are returned as local coordinates.
     *
     * \param  lattice pointer to the lattice
     * \param  ind     global coordinates of the lattice site (Input)
     * \param  grid     local coordinates of the lattice site (Output)
     * \return         index of the node for the lattice site
     */
    int map_lattice_to_node(int *ind, int *grid);

    /********************** static Functions **********************/

    /** Map a spatial position to the surrounding lattice sites.
     *
     * This function takes a global spatial position and determines the
     * surrounding elementary cell of the lattice for this position.
     * The distance fraction in each direction is also calculated.
     * <br><em>Remarks:</em>
     * <ul>
     * <li>The spatial position is given in global coordinates.</li>
     * <li>The lattice sites of the elementary cell are returned as local indices</li>
     * </ul>
     * \param pos        spatial position (Input)
     * \param ind        global index of the lower left lattice site (Output)
     * \param delta      distance fraction of pos from the surrounding
     *                   elementary cell, 6 directions (Output)
     * \param tmp_agrid  lattice mesh distance
     */
    static void map_position_to_lattice_global (double pos[3], int ind[3], double delta[6], double tmp_agrid);
};

#endif /* LATTICE_HPP */
