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
#ifndef _FFT_H
#define _FFT_H
/** \file
 *
 *  Routines, row decomposition, data structures and communication for the
 * 3D-FFT.
 *
 *  The 3D-FFT is split into 3 ond dimensional FFTs. The data is
 *  distributed in such a way, that for the actual direction of the
 *  FFT each node has a certain number of rows for which it performs a
 *  1D-FFT. After performing the FFT on that direction the data is
 *  redistributed.
 *
 *  For simplicity at the moment I have implemented a full complex to
 *  complex FFT (even though a real to complex FFT would be
 *  sufficient)
 *
 *  \todo Combine the forward and backward structures.
 *  \todo The packing routines could be moved to utils.hpp when they are needed
 * elsewhere.
 *
 *  For more information about FFT usage, see \ref fft.cpp "fft.cpp".
 */

#include "config.hpp"
#if defined(P3M) || defined(DP3M)

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>
#include <fftw3.h>

/************************************************
 * data types
 ************************************************/

/** Aligned allocator for fft data. */
template <class T> struct fft_allocator {
  typedef T value_type;
  fft_allocator() noexcept = default; // default ctor not required
  template <class U> explicit fft_allocator(const fft_allocator<U> &) {}
  template <class U> bool operator==(const fft_allocator<U> &) const {
    return true;
  }
  template <class U> bool operator!=(const fft_allocator<U> &) const {
    return false;
  }

  T *allocate(const size_t n) const {
    if (n == 0) {
      return nullptr;
    }
    if (n > static_cast<size_t>(-1) / sizeof(T)) {
      throw std::bad_array_new_length();
    }
    void *const pv = fftw_malloc(n * sizeof(T));
    if (!pv) {
      throw std::bad_alloc();
    }
    return static_cast<T *>(pv);
  }
  void deallocate(T *const p, size_t) const noexcept { fftw_free(p); }
};

template <class T> using fft_vector = std::vector<T, fft_allocator<T>>;

/** Structure for performing a 1D FFT.
 *
 *  This includes the information about the redistribution of the 3D
 *  FFT *grid before the actual FFT.
 */
struct fft_forw_plan {
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

  /** size of local mesh before communication. */
  int old_mesh[3];
  /** size of local mesh after communication, also used for actual FFT. */
  int new_mesh[3];
  /** lower left point of local FFT mesh in global FFT mesh coordinates. */
  int start[3];
  /** size of new mesh (number of mesh points). */
  int new_size;

  /** group of nodes which have to communicate with each other. */
  std::vector<int> group;

  /** packing function for send blocks. */
  void (*pack_function)(double const *const, double *const, int const *,
                        int const *, int const *, int);
  /** Send block specification. 6 integers for each node: start[3], size[3]. */
  std::vector<int> send_block;
  /** Send block communication sizes. */
  std::vector<int> send_size;
  /** Recv block specification. 6 integers for each node: start[3], size[3]. */
  std::vector<int> recv_block;
  /** Recv block communication sizes. */
  std::vector<int> recv_size;
  /** size of send block elements. */
  int element;
};

/** Additional information for backwards FFT.*/
struct fft_back_plan {
  /** plan direction. (e.g. fftw macro)*/
  int dir;
  /** plan for fft. */
  fftw_plan our_fftw_plan;

  /** packing function for send blocks. */
  void (*pack_function)(double const *const, double *const, int const *,
                        int const *, int const *, int);
};

/** Information about the three one dimensional FFTs and how the nodes
 *  have to communicate inbetween.
 *
 *  @note FFT numbering starts with 1 for technical reasons (because we have 4
 *        node grids, the index 0 is used for the real space charge assignment
 *        grid).
 */
struct fft_data_struct {
  /** Information for forward FFTs. */
  fft_forw_plan plan[4];
  /** Information for backward FFTs. */
  fft_back_plan back[4];

  /** Whether FFT is initialized or not. */
  bool init_tag = false;

  /** Maximal size of the communication buffers. */
  int max_comm_size = 0;

  /** Maximal local mesh size. */
  int max_mesh_size = 0;

  /** send buffer. */
  std::vector<double> send_buf;
  /** receive buffer. */
  std::vector<double> recv_buf;
  /** Buffer for receive data. */
  fft_vector<double> data_buf;
};

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize everything connected to the 3D-FFT.
 *
 *  \param ca_mesh_dim     Local CA mesh dimensions.
 *  \param ca_mesh_margin  Local CA mesh margins.
 *  \param global_mesh_dim Global CA mesh dimensions.
 *  \param global_mesh_off Global CA mesh offset.
 *  \param ks_pnum         Number of permutations in k-space.
 *  \param fft             FFT plan.
 *  \param grid            Number of nodes in each spatial dimension.
 *  \param comm            MPI communicator.
 *  \return Maximal size of local fft mesh (needed for allocation of ca_mesh).
 */
int fft_init(const Utils::Vector3i &ca_mesh_dim, int const *ca_mesh_margin,
             int *global_mesh_dim, double *global_mesh_off, int *ks_pnum,
             fft_data_struct &fft, const Utils::Vector3i &grid,
             const boost::mpi::communicator &comm);

/** Perform an in-place forward 3D FFT.
 *  \warning The content of \a data is overwritten.
 *  \param[in,out] data  Mesh.
 *  \param fft           FFT plan.
 *  \param comm          MPI communicator
 */
void fft_perform_forw(double *data, fft_data_struct &fft,
                      const boost::mpi::communicator &comm);

/** Perform an in-place backward 3D FFT.
 *  \warning The content of \a data is overwritten.
 *  \param[in,out] data   Mesh.
 *  \param check_complex  Throw an error if the complex component is non-zero.
 *  \param fft            FFT plan.
 *  \param comm           MPI communicator.
 */
void fft_perform_back(double *data, bool check_complex, fft_data_struct &fft,
                      const boost::mpi::communicator &comm);

/** pack a block (size[3] starting at start[3]) of an input 3d-grid
 *  with dimension dim[3] into an output 3d-block with dimension size[3].
 *
 *    The block with dimensions (size[0], size[1], size[2]) is stored
 *    in 'row-major-order' or 'C-order', that means the first index is
 *    changing slowest when running through the linear array. The
 *    element (i0 (slow), i1 (mid), i2 (fast)) has the linear index
 *    li = i2 + size[2] * (i1 + (size[1]*i0))
 *
 *  \param[in]  in      pointer to input 3d-grid.
 *  \param[out] out     pointer to output 3d-grid (block).
 *  \param[in]  start   start index of the block in the in-grid.
 *  \param[in]  size    size of the block (=dimension of the out-grid).
 *  \param[in]  dim     size of the in-grid.
 *  \param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void fft_pack_block(double const *in, double *out, int const start[3],
                    int const size[3], int const dim[3], int element);

/** unpack a 3d-grid input block (size[3]) into an output 3d-grid
 *  with dimension dim[3] at start position start[3].
 *
 *  see also \ref fft_pack_block.
 *
 *  \param[in]  in      pointer to input 3d-grid.
 *  \param[out] out     pointer to output 3d-grid (block).
 *  \param[in]  start   start index of the block in the in-grid.
 *  \param[in]  size    size of the block (=dimension of the out-grid).
 *  \param[in]  dim     size of the in-grid.
 *  \param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void fft_unpack_block(double const *in, double *out, int const start[3],
                      int const size[3], int const dim[3], int element);

/*@}*/
#endif

#endif
