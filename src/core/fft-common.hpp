/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#ifndef _FFT_COMMON_H
#define _FFT_COMMON_H

#include "config.hpp"
#if defined(P3M) || defined(DP3M)

#include <fftw3.h>

/************************************************
 * data types
 ************************************************/

/** Structure for performing a 1D FFT.  
 *
 *  This includes the information about the redistribution of the 3D
 *  FFT *grid before the actual FFT.  
*/
typedef struct {
  /** plan direction: 0 = Forward FFT, 1 = Backward FFT. */
  int dir;
  /** row direction of that FFT. */
  int row_dir;
  /** permutations from normal coordinate system. */
  int n_permute;
  /** number of 1D FFTs. */ 
  int n_ffts;
  /** plan for fft. */
  fftw_plan our_fftw_plan;
  /** function for fft. */
  void (*fft_function)(fftw_plan);

  /** size of local mesh before communication. */
  int old_mesh[3];
  /** size of local mesh after communication, also used for actual FFT. */
  int new_mesh[3];
  /** lower left point of local FFT mesh in global FFT mesh coordinates. */
  int start[3];
  /** size of new mesh (number of mesh points). */
  int new_size;

  /** number of nodes which have to communicate with each other. */ 
  int g_size;
  /** group of nodes which have to communicate with each other. */ 
  int *group;

  /** packing function for send blocks. */
  void (*pack_function)(double*, double*, int*, int*, int*, int);
  /** Send block specification. 6 integers for each node: start[3], size[3]. */ 
  int *send_block;
  /** Send block communication sizes. */ 
  int *send_size;
  /** Recv block specification. 6 integers for each node: start[3], size[3]. */ 
  int *recv_block;
  /** Recv block communication sizes. */ 
  int *recv_size;
  /** size of send block elements. */
  int element;
} fft_forw_plan;

/** Additional information for backwards FFT.*/
typedef struct {
  /** plan direction. (e.g. fftw makro)*/
  int dir;
  /** plan for fft. */
  fftw_plan our_fftw_plan;
  /** function for fft. */
  void (*fft_function)(fftw_plan);

  /** packing function for send blocks. */
  void (*pack_function)(double*, double*, int*, int*, int*, int);
} fft_back_plan;

typedef struct {
  /** Information about the three one dimensional FFTs and how the nodes
   *  have to communicate in between.
   *
   * NOTE: FFT numbering starts with 1 for technical reasons (because we
   *       have 4 node grids, the index 0 is used for the real space
   *       charge assignment grid).  */
  fft_forw_plan plan[4];
  /** Information for Back FFTs (see fft.plan). */
  fft_back_plan back[4];

  /** Whether FFT is initialized or not. */
  int init_tag;

  /** Maximal size of the communication buffers. */
  int max_comm_size;

  /** Maximal local mesh size. */
  int max_mesh_size;

  /** send buffer. */
  double *send_buf;
  /** receive buffer. */
  double *recv_buf;
  /** Buffer for receive data. */
  double *data_buf;
} fft_data_struct;

/************************************************
 * DEFINES
 ************************************************/

/* MPI tags for the fft communications: */
/** Tag for communication in fft_init() */
#define REQ_FFT_INIT   300
/** Tag for communication in forw_grid_comm() */
#define REQ_FFT_FORW   301
/** Tag for communication in back_grid_comm() */
#define REQ_FFT_BACK   302
/* Tag for wisdom file I/O */
#  define FFTW_FAILURE 0

/** Initialize FFT data structure. */
void fft_common_pre_init(fft_data_struct *fft);

/** This ugly function does the bookkepping which nodes have to
 *  communicate to each other, when you change the node grid.
 *  Changing the domain decomposition requieres communication. This
 *  function finds (hopefully) the best way to do this. As input it
 *  needs the two grids (grid1, grid2) and a linear list (node_list1)
 *  with the node identities for grid1. The linear list (node_list2)
 *  for the second grid is calculated. For the communication group of
 *  the calling node it calculates a list (group) with the node
 *  identities and the positions (pos1, pos2) of that nodes in grid1
 *  and grid2. The return value is the size of the communication
 *  group. It gives -1 if the two grids do not fit to each other
 *  (grid1 and grid2 have to be component wise multiples of each
 *  other. see e.g. \ref calc_2d_grid in \ref grid.cpp for how to do
 *  this.).
 *
 * \param grid1       The node grid you start with (Input).
 * \param grid2       The node grid you want to have (Input).
 * \param node_list1  Linear node index list for grid1 (Input).
 * \param node_list2  Linear node index list for grid2 (Output).
 * \param group       communication group (node identity list) for the calling node  (Output).
 * \param pos        positions of the nodes in in grid2 (Output).
 * \param my_pos      position of this_node in  grid2.
 * \return Size of the communication group (Output of course!).  */
int fft_find_comm_groups(int grid1[3], int grid2[3], int *node_list1, int *node_list2, 
			 int *group, int *pos, int *my_pos);


/** Calculate the local fft mesh.  Calculate the local mesh (loc_mesh)
 *  of a node at position (n_pos) in a node grid (n_grid) for a global
 *  mesh of size (mesh) and a mesh offset (mesh_off (in mesh units))
 *  and store also the first point (start) of the local mesh.
 *
 * \return size     number of mesh points in local mesh.
 * \param  n_pos    Position of the node in n_grid.
 * \param  n_grid   node grid.
 * \param  mesh     global mesh dimensions.
 * \param  mesh_off global mesh offset (see \ref p3m_data_struct).
 * \param  loc_mesh local mesh dimension (output).
 * \param  start    first point of local mesh in global mesh (output).
*/
int fft_calc_local_mesh(int n_pos[3], int n_grid[3], int mesh[3], double mesh_off[3], 
			int loc_mesh[3], int start[3]);

/** Calculate a send (or recv.) block for grid communication during a
 *  decomposition change.  Calculate the send block specification
 *  (block = lower left corner and upper right corner) which a node at
 *  position (pos1) in the actual node grid (grid1) has to send to
 *  another node at position (pos2) in the desired node grid
 *  (grid2). The global mesh, subject to communication, is specified
 *  via its size (mesh) and its mesh offset (mesh_off (in mesh
 *  units)).
 *
 *  For the calculation of a receive block you have to change the arguments in the following way: <br>
 *  pos1  - position of receiving node in the desired node grid. <br>
 *  grid1 - desired node grid. <br>
 *  pos2  - position of the node you intend to receive the data from in the actual node grid. <br>
 *  grid2 - actual node grid.  <br>
 *
 *  \return          size of the send block.
 *  \param  pos1     Position of send node in grid1.
 *  \param  grid1    node grid 1.
 *  \param  pos2     Position of recv node in grid2.
 *  \param  grid2    node grid 2.
 *  \param  mesh     global mesh dimensions.
 *  \param  mesh_off global mesh offset (see \ref p3m_data_struct).
 *  \param  block    send block specification.
*/
int fft_calc_send_block(int pos1[3], int grid1[3], int pos2[3], int grid2[3], 
			int mesh[3], double mesh_off[3], int block[6]);


/** pack a block (size[3] starting at start[3]) of an input 3d-grid
 *  with dimension dim[3] into an output 3d-block with dimension size[3].
 *
 *    The block with dimensions (size[0], size[1], size[2]) is stored
 *    in 'row-major-order' or 'C-order', that means the first index is
 *    changing slowest when running through the linear array. The
 *    element (i0 (slow), i1 (mid), i2 (fast)) has the linear index 
 *    li = i2 + size[2] * (i1 + (size[1]*i0)) 
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void fft_pack_block(double *in, double *out, int start[3], int size[3], 
		    int dim[3], int element);

/** pack a block with dimensions (size[0] * size[1] * aize[2]) starting
 *  at start[3] of an input 3d-grid with dimension dim[3] into an
 *  output 3d-grid with dimensions (size[2] * size[0] * size[1]) with
 *  a simulatanous one-fold permutation of the indices.
 *
 * The permutation is defined as: 
 * slow_in -> fast_out, mid_in ->slow_out, fast_in -> mid_out
 *
 * An element (i0_in , i1_in , i2_in ) is then 
 * (i0_out = i1_in-start[1], i1_out = i2_in-start[2], i2_out = i0_in-start[0]) and
 * for the linear indices we have:                              <br>
 * li_in = i2_in + size[2] * (i1_in + (size[1]*i0_in))          <br>
 * li_out = i2_out + size[0] * (i1_out + (size[2]*i0_out)) 
 *
 * For index definition see \ref fft_pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void fft_pack_block_permute1(double *in, double *out, int start[3], int size[3], 
			     int dim[3], int element);

/** pack a block with dimensions (size[0] * size[1] * aize[2]) starting
 *  at start[3] of an input 3d-grid with dimension dim[3] into an
 *  output 3d-grid with dimensions (size[2] * size[0] * size[1]), this
 *  is a simulatanous two-fold permutation of the indices.
 *
 * The permutation is defined as: 
 * slow_in -> mid_out, mid_in ->fast_out, fast_in -> slow_out
 *
 * An element (i0_in , i1_in , i2_in ) is then 
 * (i0_out = i2_in-start[2], i1_out = i0_in-start[0], i2_out = i1_in-start[1]) and
 * for the linear indices we have:                              <br>
 * li_in = i2_in + size[2] * (i1_in + (size[1]*i0_in))          <br>
 * li_out = i2_out + size[0] * (i1_out + (size[2]*i0_out)) 
 *
 * For index definition see \ref fft_pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void fft_pack_block_permute2(double *in, double *out, int start[3], int size[3], 
			 int dim[3],int element);


/** unpack a 3d-grid input block (size[3]) into an output 3d-grid
 *  with dimension dim[3] at start position start[3].
 *
 *  see also \ref fft_pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void fft_unpack_block(double *in, double *out, int start[3], int size[3], 
		      int dim[3], int element);

/** Debug function to print global fft mesh. 
    Print a globaly distributed mesh contained in data. Element size is element. 
 * \param plan     fft/communication plan (see \ref fft_forw_plan).
 * \param data     mesh data.
 * \param element  element size.
 * \param num      element index to print.
*/
void fft_print_global_fft_mesh(fft_forw_plan plan, double *data, int element, int num);

/** Debug function to print fft_forw_plan structure. 
 * \param pl fft/communication plan (see \ref fft_forw_plan).
 */
void fft_print_fft_plan(fft_forw_plan pl);


#endif /* defined(P3M) || defined(DP3M) */
#endif /* _FFT_COMMON_H */
