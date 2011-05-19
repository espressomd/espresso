/*
  Copyright (C) 2010,2011 The ESPResSo project
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
#ifndef _FFT_H
#define _FFT_H
/** \file fft.h
 *
 *  Routines, row decomposition, data structures and communication for the 3D-FFT. 
 *
 *  The 3D-FFT is split into 3 ond dimensional FFTs. The data is
 *  distributed in such a way, that for the actual direction of the
 *  FFT each node has a certain number of rows for which it performs a
 *  1D-FFT. After performing the FFT on theat direction the data is
 *  redistributed.
 *
 *  For simplicity at the moment I have implemented a full complex to
 *  complex FFT (even though a real to complex FFT would be
 *  sufficient)
 *
 *  \todo Combine the forward and backward structures.
 *  \todo The packing routines could be moved to utils.h when they are needed elsewhere.
 *
 *  For more information about FFT usage, see \ref fft.c "fft.c".  
*/

#include "utils.h"

#if defined(P3M) || defined(DP3M)

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
  void *fft_plan;
  /** function for fft. */
  void (*fft_function)();

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
  void (*pack_function)();
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
  void *fft_plan;
  /** function for fft. */
  void (*fft_function)();

  /** packing function for send blocks. */
  void (*pack_function)(); 
} fft_back_plan;


/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Information about the three one dimensional FFTs and how the nodes
 *  have to communicate in between.
 *
 * NOTE: FFT numbering starts with 1 for technical reasons (because we
 *       have 4 node grids, the index 0 is used for the real space
 *       charge assignment grid).  */

#ifdef P3M
extern fft_forw_plan fft_plan[4];
#endif

#ifdef DP3M
extern fft_forw_plan Dfft_plan[4];
#endif


/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize some arrays connected to the 3D-FFT. */
void  fft_pre_init();

#ifdef P3M
/** Initialize everything connected to the 3D-FFT.

 * \return Maximal size of local fft mesh (needed for allocation of ca_mesh).
 * \param data           Pointer Pounter to data array.
 * \param ca_mesh_dim    Pointer to CA mesh dimensions.
 * \param ca_mesh_margin Pointer to CA mesh margins.
 * \param ks_pnum        Pointer to number of permutations in k-space.
 */
int fft_init(double **data, int *ca_mesh_dim, int *ca_mesh_margin, int *ks_pnum);

/** perform the forward 3D FFT.
    The assigned charges are in \a data. The result is also stored in \a data.
    \warning The content of \a data is overwritten.
    \param data Mesh.
*/
void fft_perform_forw(double *data);

/** perform the backward 3D FFT.
    \warning The content of \a data is overwritten.
    \param data Mesh.
*/
void fft_perform_back(double *data);

#endif

#ifdef DP3M
/** Initialize everything connected to the 3D-FFT related to the dipole-dipole.

 * \return Maximal size of local fft mesh (needed for allocation of ca_mesh).
 * \param data           Pointer Pounter to data array.
 * \param ca_mesh_dim    Pointer to CA mesh dimensions.
 * \param ca_mesh_margin Pointer to CA mesh margins.
 * \param ks_pnum        Pointer to number of permutations in k-space.
 */
int Dfft_init(double **data, int *ca_mesh_dim, int *ca_mesh_margin, int *ks_pnum);

/** perform the forward 3D FFT for meshes related to the magnetic dipole-dipole interaction.
    The assigned charges are in \a data. The result is also stored in \a data.
    \warning The content of \a data is overwritten.
    \param data DMesh.
*/
void Dfft_perform_forw(double *data);

/** perform the backward 3D FFT for meshes related to the magnetic dipole-dipole interaction.
    \warning The content of \a data is overwritten.
    \param data DMesh.
*/
void Dfft_perform_back(double *data);

#endif






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
void pack_block(double *in, double *out, int start[3], int size[3], 
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
 * For index definition see \ref pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void pack_block_permute1(double *in, double *out, int start[3], int size[3], 
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
 * For index definition see \ref pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void pack_block_permute2(double *in, double *out, int start[3], int size[3], 
			 int dim[3],int element);


/** unpack a 3d-grid input block (size[3]) into an output 3d-grid
 *  with dimension dim[3] at start position start[3].
 *
 *  see also \ref pack_block.
 *
 *  \param in      pointer to input 3d-grid.
 *  \param out     pointer to output 3d-grid (block).
 *  \param start   start index of the block in the in-grid.
 *  \param size    size of the block (=dimension of the out-grid).
 *  \param dim     size of the in-grid.
 *  \param element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void unpack_block(double *in, double *out, int start[3], int size[3], 
		  int dim[3], int element);

/** Debug function to print global fft mesh. 
    Print a globaly distributed mesh contained in data. Element size is element. 
 * \param plan     fft/communication plan (see \ref fft_forw_plan).
 * \param data     mesh data.
 * \param element  element size.
 * \param num      element index to print.
*/
void print_global_fft_mesh(fft_forw_plan plan, double *data, int element, int num);

/*@}*/
#endif

#endif
