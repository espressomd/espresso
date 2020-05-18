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
/** \file
 *
 *  Routines, row decomposition, data structures and communication for the
 *  3D-FFT.
 *
 */

#include "fft.hpp"

#if defined(P3M) || defined(DP3M)

#include <utils/math/permute_ifield.hpp>
using Utils::permute_ifield;
#include <utils/index.hpp>
using Utils::get_linear_index;

#include <fftw3.h>
#include <mpi.h>

#include <cstdio>
#include <cstring>

#include <utils/Span.hpp>

/************************************************
 * DEFINES
 ************************************************/

/** @name MPI tags for FFT communication */
/*@{*/
/** Tag for communication in fft_init() */
/** Tag for communication in forw_grid_comm() */
#define REQ_FFT_FORW 301
/** Tag for communication in back_grid_comm() */
#define REQ_FFT_BACK 302
/*@}*/

namespace {
/** This ugly function does the bookkeeping: which nodes have to
 *  communicate to each other, when you change the node grid.
 *  Changing the domain decomposition requires communication. This
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
 *  \param[in]  grid1       The node grid you start with.
 *  \param[in]  grid2       The node grid you want to have.
 *  \param[in]  node_list1  Linear node index list for grid1.
 *  \param[out] node_list2  Linear node index list for grid2.
 *  \param[out] pos         Positions of the nodes in grid2
 *  \param[out] my_pos      Position of comm.rank() in grid2.
 *  \return Size of the communication group.
 */
boost::optional<std::vector<int>>
find_comm_groups(Utils::Vector3i const &grid1, Utils::Vector3i const &grid2,
                 Utils::Span<const int> node_list1, Utils::Span<int> node_list2,
                 Utils::Span<int> pos, Utils::Span<int> my_pos,
                 boost::mpi::communicator const &comm) {
  int i;
  /* communication group cell size on grid1 and grid2 */
  int s1[3], s2[3];
  /* The communication group cells build the same super grid on grid1 and grid2
   */
  int ds[3];
  /* communication group size */
  int g_size = 1;
  /* comm. group cell index */
  int gi[3];
  /* position of a node in a grid */
  int p1[3], p2[3];
  /* node identity */
  int n;
  /* comm.rank() position in the communication group. */
  int c_pos = -1;
  /* flag for group identification */
  int my_group = 0;

  /* calculate dimension of comm. group cells for both grids */
  if ((grid1[0] * grid1[1] * grid1[2]) != (grid2[0] * grid2[1] * grid2[2]))
    return boost::none; /* unlike number of nodes */
  for (i = 0; i < 3; i++) {
    s1[i] = grid1[i] / grid2[i];
    if (s1[i] == 0)
      s1[i] = 1;
    else if (grid1[i] != grid2[i] * s1[i])
      return boost::none; /* grids do not match!!! */

    s2[i] = grid2[i] / grid1[i];
    if (s2[i] == 0)
      s2[i] = 1;
    else if (grid2[i] != grid1[i] * s2[i])
      return boost::none; /* grids do not match!!! */

    ds[i] = grid2[i] / s2[i];
    g_size *= s2[i];
  }

  std::vector<int> group(g_size);

  /* calc node_list2 */
  /* loop through all comm. group cells */
  for (gi[2] = 0; gi[2] < ds[2]; gi[2]++)
    for (gi[1] = 0; gi[1] < ds[1]; gi[1]++)
      for (gi[0] = 0; gi[0] < ds[0]; gi[0]++) {
        /* loop through all nodes in that comm. group cell */
        for (i = 0; i < g_size; i++) {
          p1[0] = (gi[0] * s1[0]) + (i % s1[0]);
          p1[1] = (gi[1] * s1[1]) + ((i / s1[0]) % s1[1]);
          p1[2] = (gi[2] * s1[2]) + (i / (s1[0] * s1[1]));

          p2[0] = (gi[0] * s2[0]) + (i % s2[0]);
          p2[1] = (gi[1] * s2[1]) + ((i / s2[0]) % s2[1]);
          p2[2] = (gi[2] * s2[2]) + (i / (s2[0] * s2[1]));

          n = node_list1[get_linear_index(p1[0], p1[1], p1[2], grid1)];
          node_list2[get_linear_index(p2[0], p2[1], p2[2], grid2)] = n;

          pos[3 * n + 0] = p2[0];
          pos[3 * n + 1] = p2[1];
          pos[3 * n + 2] = p2[2];
          if (my_group == 1)
            group[i] = n;
          if (n == comm.rank() && my_group == 0) {
            my_group = 1;
            c_pos = i;
            my_pos[0] = p2[0];
            my_pos[1] = p2[1];
            my_pos[2] = p2[2];
            i = -1; /* restart the loop */
          }
        }
        my_group = 0;
      }

  /* permute comm. group according to the nodes position in the group */
  /* This is necessary to have matching node pairs during communication! */
  while (c_pos > 0) {
    n = group[g_size - 1];
    for (i = g_size - 1; i > 0; i--)
      group[i] = group[i - 1];
    group[0] = n;
    c_pos--;
  }
  return group;
}

/** Calculate the local fft mesh. Calculate the local mesh (@p loc_mesh)
 *  of a node at position (@p n_pos) in a node grid (@p n_grid) for a global
 *  mesh of size (@p mesh) and a mesh offset (@p mesh_off (in mesh units))
 *  and store also the first point (@p start) of the local mesh.
 *
 * \param[in]  n_pos    Position of the node in @p n_grid.
 * \param[in]  n_grid   node grid.
 * \param[in]  mesh     global mesh dimensions.
 * \param[in]  mesh_off global mesh offset (see \ref p3m_data_struct).
 * \param[out] loc_mesh local mesh dimension.
 * \param[out] start    first point of local mesh in global mesh.
 * \return Number of mesh points in local mesh.
 */
int calc_local_mesh(const int *n_pos, const int *n_grid, const int *mesh,
                    const double *mesh_off, int *loc_mesh, int *start) {
  int last[3], size = 1;

  for (int i = 0; i < 3; i++) {
    start[i] =
        (int)ceil((mesh[i] / (double)n_grid[i]) * n_pos[i] - mesh_off[i]);
    last[i] = (int)floor((mesh[i] / (double)n_grid[i]) * (n_pos[i] + 1) -
                         mesh_off[i]);
    /* correct round off errors */
    if ((mesh[i] / (double)n_grid[i]) * (n_pos[i] + 1) - mesh_off[i] - last[i] <
        1.0e-15)
      last[i]--;
    if (1.0 + (mesh[i] / (double)n_grid[i]) * n_pos[i] - mesh_off[i] -
            start[i] <
        1.0e-15)
      start[i]--;
    loc_mesh[i] = last[i] - start[i] + 1;
    size *= loc_mesh[i];
  }
  return size;
}

/** Calculate a send (or recv.) block for grid communication during a
 *  decomposition change. Calculate the send block specification
 *  (block = lower left corner and upper right corner) which a node at
 *  position (@p pos1) in the actual node grid (@p grid1) has to send to
 *  another node at position (@p pos2) in the desired node grid (@p grid2).
 *  The global mesh, subject to communication, is specified via its size
 *  (@p mesh) and its mesh offset (@p mesh_off (in mesh units)).
 *
 *  For the calculation of a receive block you have to change the arguments in
 *  the following way:
 *  - @p pos1: position of receiving node in the desired node grid.
 *  - @p grid1: desired node grid.
 *  - @p pos2: position of the node you intend to receive the data from in the
 *    actual node grid.
 *  - @p grid2: actual node grid.
 *
 *  \param[in]  pos1     Position of send node in @p grid1.
 *  \param[in]  grid1    node grid 1.
 *  \param[in]  pos2     Position of recv node in @p grid2.
 *  \param[in]  grid2    node grid 2.
 *  \param[in]  mesh     global mesh dimensions.
 *  \param[in]  mesh_off global mesh offset (see \ref p3m_data_struct).
 *  \param[out] block    send block specification.
 *  \return Size of the send block.
 */
int calc_send_block(const int *pos1, const int *grid1, const int *pos2,
                    const int *grid2, const int *mesh, const double *mesh_off,
                    int *block) {
  int size = 1;
  int mesh1[3], first1[3], last1[3];
  int mesh2[3], first2[3], last2[3];

  calc_local_mesh(pos1, grid1, mesh, mesh_off, mesh1, first1);
  calc_local_mesh(pos2, grid2, mesh, mesh_off, mesh2, first2);

  for (int i = 0; i < 3; i++) {
    last1[i] = first1[i] + mesh1[i] - 1;
    last2[i] = first2[i] + mesh2[i] - 1;
    block[i] = std::max(first1[i], first2[i]) - first1[i];
    block[i + 3] = (std::min(last1[i], last2[i]) - first1[i]) - block[i] + 1;
    size *= block[i + 3];
  }
  return size;
}

/** Pack a block with dimensions <tt>size[0] * size[1] * size[2]</tt> starting
 *  at @p start of an input 3D-grid with dimension @p dim into an output
 *  3D-grid with dimensions <tt>size[2] * size[0] * size[1]</tt> with
 *  a simultaneous one-fold permutation of the indices. The permutation is
 *  defined as: slow_in -> fast_out, mid_in ->slow_out, fast_in -> mid_out.
 *
 *  An element <tt>(i0_in, i1_in, i2_in)</tt> is then
 *  <tt>(i0_out = i1_in-start[1], i1_out = i2_in-start[2],
 *  i2_out = i0_in-start[0])</tt> and for the linear indices we have:
 *  - <tt>li_in = i2_in + size[2] * (i1_in + (size[1]*i0_in))</tt>
 *  - <tt>li_out = i2_out + size[0] * (i1_out + (size[2]*i0_out))</tt>
 *
 *  For index definition see \ref fft_pack_block.
 *
 *  \param[in]  in      input 3D-grid.
 *  \param[out] out     output 3D-grid (block).
 *  \param[in]  start   start index of the block in the in-grid.
 *  \param[in]  size    size of the block (=dimension of the out-grid).
 *  \param[in]  dim     size of the in-grid.
 *  \param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void pack_block_permute1(double const *const in, double *const out,
                         const int *start, const int *size, const int *dim,
                         int element) {
  /* slow, mid and fast changing indices for input grid */
  int s, m, f, e;
  /* linear index of in grid, linear index of out grid */
  int li_in, li_out = 0;
  /* offsets for indices in input grid */
  int m_in_offset, s_in_offset;
  /* offset for mid changing indices of output grid */
  int m_out_offset;

  m_in_offset = element * (dim[2] - size[2]);
  s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  m_out_offset = (element * size[0]) - element;
  li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (s = 0; s < size[0]; s++) { /* fast changing out */
    li_out = element * s;
    for (m = 0; m < size[1]; m++) {   /* slow changing out */
      for (f = 0; f < size[2]; f++) { /* mid  changing out */
        for (e = 0; e < element; e++)
          out[li_out++] = in[li_in++];
        li_out += m_out_offset;
      }
      li_in += m_in_offset;
    }
    li_in += s_in_offset;
  }
}

/** Pack a block with dimensions <tt>size[0] * size[1] * size[2]</tt> starting
 *  at @p start of an input 3D-grid with dimension @p dim into an output
 *  3D-grid with dimensions <tt>size[2] * size[0] * size[1]</tt> with
 *  a simultaneous two-fold permutation of the indices. The permutation is
 *  defined as: slow_in -> mid_out, mid_in ->fast_out, fast_in -> slow_out.
 *
 *  An element <tt>(i0_in, i1_in, i2_in)</tt> is then
 *  <tt>(i0_out = i2_in-start[2], i1_out = i0_in-start[0],
 *  i2_out = i1_in-start[1])</tt> and for the linear indices we have:
 *  - <tt>li_in = i2_in + size[2] * (i1_in + (size[1]*i0_in))</tt>
 *  - <tt>li_out = i2_out + size[0] * (i1_out + (size[2]*i0_out))</tt>
 *
 *  For index definition see \ref fft_pack_block.
 *
 *  \param[in]  in      input 3D-grid.
 *  \param[out] out     output 3D-grid (block).
 *  \param[in]  start   start index of the block in the in-grid.
 *  \param[in]  size    size of the block (=dimension of the out-grid).
 *  \param[in]  dim     size of the in-grid.
 *  \param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
void pack_block_permute2(double const *const in, double *const out,
                         const int *start, const int *size, const int *dim,
                         int element) {
  /* slow, mid and fast changing indices for input grid */
  int s, m, f, e;
  /* linear index of in grid, linear index of out grid */
  int li_in, li_out = 0;
  /* offsets for indices in input grid */
  int m_in_offset, s_in_offset;
  /* offset for slow changing index of output grid */
  int s_out_offset;
  /* start index for mid changing index of output grid */
  int m_out_start;

  m_in_offset = element * (dim[2] - size[2]);
  s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  s_out_offset = (element * size[0] * size[1]) - element;
  li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (s = 0; s < size[0]; s++) { /* mid changing out */
    m_out_start = element * (s * size[1]);
    for (m = 0; m < size[1]; m++) { /* fast changing out */
      li_out = m_out_start + element * m;
      for (f = 0; f < size[2]; f++) { /* slow  changing out */
        for (e = 0; e < element; e++)
          out[li_out++] = in[li_in++];
        li_out += s_out_offset;
      }
      li_in += m_in_offset;
    }
    li_in += s_in_offset;
  }
}

/** Communicate the grid data according to the given forward FFT plan.
 *  \param plan   FFT communication plan.
 *  \param in     input mesh.
 *  \param out    output mesh.
 *  \param fft    FFT communication plan.
 *  \param comm   MPI communicator.
 */
void forw_grid_comm(fft_forw_plan plan, const double *in, double *out,
                    fft_data_struct &fft,
                    const boost::mpi::communicator &comm) {
  for (int i = 0; i < plan.group.size(); i++) {
    plan.pack_function(in, fft.send_buf.data(), &(plan.send_block[6 * i]),
                       &(plan.send_block[6 * i + 3]), plan.old_mesh,
                       plan.element);

    if (plan.group[i] != comm.rank()) {
      MPI_Sendrecv(fft.send_buf.data(), plan.send_size[i], MPI_DOUBLE,
                   plan.group[i], REQ_FFT_FORW, fft.recv_buf.data(),
                   plan.recv_size[i], MPI_DOUBLE, plan.group[i], REQ_FFT_FORW,
                   comm, MPI_STATUS_IGNORE);
    } else { /* Self communication... */
      std::swap(fft.send_buf, fft.recv_buf);
    }
    fft_unpack_block(fft.recv_buf.data(), out, &(plan.recv_block[6 * i]),
                     &(plan.recv_block[6 * i + 3]), plan.new_mesh,
                     plan.element);
  }
}

/** Communicate the grid data according to the given backward FFT plan.
 *  \param plan_f Forward FFT plan.
 *  \param plan_b Backward FFT plan.
 *  \param in     input mesh.
 *  \param out    output mesh.
 *  \param fft    FFT communication plan.
 *  \param comm   MPI communicator.
 */
void back_grid_comm(fft_forw_plan plan_f, fft_back_plan plan_b,
                    const double *in, double *out, fft_data_struct &fft,
                    const boost::mpi::communicator &comm) {
  /* Back means: Use the send/receive stuff from the forward plan but
     replace the receive blocks by the send blocks and vice
     versa. Attention then also new_mesh and old_mesh are exchanged */

  for (int i = 0; i < plan_f.group.size(); i++) {
    plan_b.pack_function(in, fft.send_buf.data(), &(plan_f.recv_block[6 * i]),
                         &(plan_f.recv_block[6 * i + 3]), plan_f.new_mesh,
                         plan_f.element);

    if (plan_f.group[i] != comm.rank()) { /* send first, receive second */
      MPI_Sendrecv(fft.send_buf.data(), plan_f.recv_size[i], MPI_DOUBLE,
                   plan_f.group[i], REQ_FFT_BACK, fft.recv_buf.data(),
                   plan_f.send_size[i], MPI_DOUBLE, plan_f.group[i],
                   REQ_FFT_BACK, comm, MPI_STATUS_IGNORE);
    } else { /* Self communication... */
      std::swap(fft.send_buf, fft.recv_buf);
    }
    fft_unpack_block(fft.recv_buf.data(), out, &(plan_f.send_block[6 * i]),
                     &(plan_f.send_block[6 * i + 3]), plan_f.old_mesh,
                     plan_f.element);
  }
}

/** Calculate 'best' mapping between a 2D and 3D grid.
 *  Required for the communication from 3D domain decomposition
 *  to 2D row decomposition.
 *  The dimensions of the 2D grid are resorted, if necessary, in a way
 *  that they are multiples of the 3D grid dimensions.
 *  \param g3d      3D grid.
 *  \param g2d      2D grid.
 *  \param mult     factors between 3D and 2D grid dimensions
 *  \return         index of the row direction [0,1,2].
 */
int map_3don2d_grid(int const g3d[3], int g2d[3], int mult[3]) {
  int row_dir = -1;
  /* trivial case */
  if (g3d[2] == 1) {
    for (int i = 0; i < 3; i++)
      mult[i] = 1;
    return 2;
  }
  if (g2d[0] % g3d[0] == 0) {
    if (g2d[1] % g3d[1] == 0) {
      row_dir = 2;
    } else if (g2d[1] % g3d[2] == 0) {
      row_dir = 1;
      g2d[2] = g2d[1];
      g2d[1] = 1;
    }
  } else if (g2d[0] % g3d[1] == 0) {
    if (g2d[1] % g3d[0] == 0) {
      row_dir = 2;
      int const tmp = g2d[0];
      g2d[0] = g2d[1];
      g2d[1] = tmp;
    } else if (g2d[1] % g3d[2] == 0) {
      row_dir = 0;
      g2d[2] = g2d[1];
      g2d[1] = g2d[0];
      g2d[0] = 1;
    }
  } else if (g2d[0] % g3d[2] == 0) {
    if (g2d[1] % g3d[0] == 0) {
      row_dir = 1;
      g2d[2] = g2d[0];
      g2d[0] = g2d[1];
      g2d[1] = 1;
    } else if (g2d[1] % g3d[1] == 0) {
      row_dir = 0;
      g2d[2] = g2d[0];
      g2d[0] = 1;
    }
  }
  for (int i = 0; i < 3; i++)
    mult[i] = g2d[i] / g3d[i];
  return row_dir;
}

/** Calculate most square 2D grid. */
void calc_2d_grid(int n, int grid[3]) {
  for (auto i = static_cast<int>(std::sqrt(n)); i >= 1; i--) {
    if (n % i == 0) {
      grid[0] = n / i;
      grid[1] = i;
      grid[2] = 1;
      return;
    }
  }
}
} // namespace

int fft_init(const Utils::Vector3i &ca_mesh_dim, int const *ca_mesh_margin,
             int *global_mesh_dim, double *global_mesh_off, int *ks_pnum,
             fft_data_struct &fft, const Utils::Vector3i &grid,
             const boost::mpi::communicator &comm) {
  int i, j;
  /* helpers */
  int mult[3];

  int n_grid[4][3];         /* The four node grids. */
  int my_pos[4][3];         /* The position of comm.rank() in the node grids. */
  std::vector<int> n_id[4]; /* linear node identity lists for the node grids. */
  std::vector<int> n_pos[4]; /* positions of nodes in the node grids. */

  int node_pos[3];
  MPI_Cart_coords(comm, comm.rank(), 3, node_pos);

  fft.max_comm_size = 0;
  fft.max_mesh_size = 0;
  for (i = 0; i < 4; i++) {
    n_id[i].resize(1 * comm.size());
    n_pos[i].resize(3 * comm.size());
  }

  /* === node grids === */
  /* real space node grid (n_grid[0]) */
  for (i = 0; i < 3; i++) {
    n_grid[0][i] = grid[i];
    my_pos[0][i] = node_pos[i];
  }
  for (i = 0; i < comm.size(); i++) {
    MPI_Cart_coords(comm, i, 3, &(n_pos[0][3 * i + 0]));
    auto const lin_ind = get_linear_index(
        n_pos[0][3 * i + 0], n_pos[0][3 * i + 1], n_pos[0][3 * i + 2],
        {n_grid[0][0], n_grid[0][1], n_grid[0][2]});
    n_id[0][lin_ind] = i;
  }

  /* FFT node grids (n_grid[1 - 3]) */
  calc_2d_grid(comm.size(), n_grid[1]);
  /* resort n_grid[1] dimensions if necessary */
  fft.plan[1].row_dir = map_3don2d_grid(n_grid[0], n_grid[1], mult);
  fft.plan[0].n_permute = 0;
  for (i = 1; i < 4; i++)
    fft.plan[i].n_permute = (fft.plan[1].row_dir + i) % 3;
  for (i = 0; i < 3; i++) {
    n_grid[2][i] = n_grid[1][(i + 1) % 3];
    n_grid[3][i] = n_grid[1][(i + 2) % 3];
  }
  fft.plan[2].row_dir = (fft.plan[1].row_dir - 1) % 3;
  fft.plan[3].row_dir = (fft.plan[1].row_dir - 2) % 3;

  /* === communication groups === */
  /* copy local mesh off real space charge assignment grid */
  for (i = 0; i < 3; i++)
    fft.plan[0].new_mesh[i] = ca_mesh_dim[i];

  for (i = 1; i < 4; i++) {
    using Utils::make_span;
    auto group = find_comm_groups(
        {n_grid[i - 1][0], n_grid[i - 1][1], n_grid[i - 1][2]},
        {n_grid[i][0], n_grid[i][1], n_grid[i][2]}, n_id[i - 1],
        make_span(n_id[i]), make_span(n_pos[i]), my_pos[i], comm);
    if (not group) {
      /* try permutation */
      std::swap(n_grid[i][(fft.plan[i].row_dir + 1) % 3],
                n_grid[i][(fft.plan[i].row_dir + 2) % 3]);

      group = find_comm_groups(
          {n_grid[i - 1][0], n_grid[i - 1][1], n_grid[i - 1][2]},
          {n_grid[i][0], n_grid[i][1], n_grid[i][2]}, make_span(n_id[i - 1]),
          make_span(n_id[i]), make_span(n_pos[i]), my_pos[i], comm);

      if (not group) {
        throw std::runtime_error("INTERNAL ERROR: fft_find_comm_groups error");
      }
    }

    fft.plan[i].group = *group;

    fft.plan[i].send_block.resize(6 * fft.plan[i].group.size());
    fft.plan[i].send_size.resize(fft.plan[i].group.size());
    fft.plan[i].recv_block.resize(6 * fft.plan[i].group.size());
    fft.plan[i].recv_size.resize(fft.plan[i].group.size());

    fft.plan[i].new_size =
        calc_local_mesh(my_pos[i], n_grid[i], global_mesh_dim, global_mesh_off,
                        fft.plan[i].new_mesh, fft.plan[i].start);
    permute_ifield(fft.plan[i].new_mesh, 3, -(fft.plan[i].n_permute));
    permute_ifield(fft.plan[i].start, 3, -(fft.plan[i].n_permute));
    fft.plan[i].n_ffts = fft.plan[i].new_mesh[0] * fft.plan[i].new_mesh[1];

    /* === send/recv block specifications === */
    for (j = 0; j < fft.plan[i].group.size(); j++) {
      /* send block: comm.rank() to comm-group-node i (identity: node) */
      int node = fft.plan[i].group[j];
      fft.plan[i].send_size[j] = calc_send_block(
          my_pos[i - 1], n_grid[i - 1], &(n_pos[i][3 * node]), n_grid[i],
          global_mesh_dim, global_mesh_off, &(fft.plan[i].send_block[6 * j]));
      permute_ifield(&(fft.plan[i].send_block[6 * j]), 3,
                     -(fft.plan[i - 1].n_permute));
      permute_ifield(&(fft.plan[i].send_block[6 * j + 3]), 3,
                     -(fft.plan[i - 1].n_permute));
      if (fft.plan[i].send_size[j] > fft.max_comm_size)
        fft.max_comm_size = fft.plan[i].send_size[j];
      /* First plan send blocks have to be adjusted, since the CA grid
         may have an additional margin outside the actual domain of the
         node */
      if (i == 1) {
        for (int k = 0; k < 3; k++)
          fft.plan[1].send_block[6 * j + k] += ca_mesh_margin[2 * k];
      }
      /* recv block: comm.rank() from comm-group-node i (identity: node) */
      fft.plan[i].recv_size[j] = calc_send_block(
          my_pos[i], n_grid[i], &(n_pos[i - 1][3 * node]), n_grid[i - 1],
          global_mesh_dim, global_mesh_off, &(fft.plan[i].recv_block[6 * j]));
      permute_ifield(&(fft.plan[i].recv_block[6 * j]), 3,
                     -(fft.plan[i].n_permute));
      permute_ifield(&(fft.plan[i].recv_block[6 * j + 3]), 3,
                     -(fft.plan[i].n_permute));
      if (fft.plan[i].recv_size[j] > fft.max_comm_size)
        fft.max_comm_size = fft.plan[i].recv_size[j];
    }

    for (j = 0; j < 3; j++)
      fft.plan[i].old_mesh[j] = fft.plan[i - 1].new_mesh[j];
    if (i == 1)
      fft.plan[i].element = 1;
    else {
      fft.plan[i].element = 2;
      for (j = 0; j < fft.plan[i].group.size(); j++) {
        fft.plan[i].send_size[j] *= 2;
        fft.plan[i].recv_size[j] *= 2;
      }
    }
  }

  /* Factor 2 for complex fields */
  fft.max_comm_size *= 2;
  fft.max_mesh_size = (ca_mesh_dim[0] * ca_mesh_dim[1] * ca_mesh_dim[2]);
  for (i = 1; i < 4; i++)
    if (2 * fft.plan[i].new_size > fft.max_mesh_size)
      fft.max_mesh_size = 2 * fft.plan[i].new_size;

  /* === pack function === */
  for (i = 1; i < 4; i++) {
    fft.plan[i].pack_function = pack_block_permute2;
  }
  (*ks_pnum) = 6;
  if (fft.plan[1].row_dir == 2) {
    fft.plan[1].pack_function = fft_pack_block;
    (*ks_pnum) = 4;
  } else if (fft.plan[1].row_dir == 1) {
    fft.plan[1].pack_function = pack_block_permute1;
    (*ks_pnum) = 5;
  }

  /* Factor 2 for complex numbers */
  fft.send_buf.resize(fft.max_comm_size);
  fft.recv_buf.resize(fft.max_comm_size);
  fft.data_buf.resize(fft.max_mesh_size);
  auto *c_data = (fftw_complex *)(fft.data_buf.data());

  /* === FFT Routines (Using FFTW / RFFTW package)=== */
  for (i = 1; i < 4; i++) {
    fft.plan[i].dir = FFTW_FORWARD;
    /* FFT plan creation.*/

    if (fft.init_tag)
      fftw_destroy_plan(fft.plan[i].our_fftw_plan);
    fft.plan[i].our_fftw_plan = fftw_plan_many_dft(
        1, &fft.plan[i].new_mesh[2], fft.plan[i].n_ffts, c_data, nullptr, 1,
        fft.plan[i].new_mesh[2], c_data, nullptr, 1, fft.plan[i].new_mesh[2],
        fft.plan[i].dir, FFTW_PATIENT);
  }

  /* === The BACK Direction === */
  /* this is needed because slightly different functions are used */
  for (i = 1; i < 4; i++) {
    fft.back[i].dir = FFTW_BACKWARD;

    if (fft.init_tag)
      fftw_destroy_plan(fft.back[i].our_fftw_plan);
    fft.back[i].our_fftw_plan = fftw_plan_many_dft(
        1, &fft.plan[i].new_mesh[2], fft.plan[i].n_ffts, c_data, nullptr, 1,
        fft.plan[i].new_mesh[2], c_data, nullptr, 1, fft.plan[i].new_mesh[2],
        fft.back[i].dir, FFTW_PATIENT);

    fft.back[i].pack_function = pack_block_permute1;
  }
  if (fft.plan[1].row_dir == 2) {
    fft.back[1].pack_function = fft_pack_block;
  } else if (fft.plan[1].row_dir == 1) {
    fft.back[1].pack_function = pack_block_permute2;
  }

  fft.init_tag = true;

  return fft.max_mesh_size;
}

void fft_perform_forw(double *data, fft_data_struct &fft,
                      const boost::mpi::communicator &comm) {
  /* ===== first direction  ===== */

  auto *c_data = (fftw_complex *)data;
  auto *c_data_buf = (fftw_complex *)fft.data_buf.data();

  /* communication to current dir row format (in is data) */
  forw_grid_comm(fft.plan[1], data, fft.data_buf.data(), fft, comm);

  /* complexify the real data array (in is fft.data_buf) */
  for (int i = 0; i < fft.plan[1].new_size; i++) {
    data[2 * i + 0] = fft.data_buf[i]; /* real value */
    data[2 * i + 1] = 0;               /* complex value */
  }
  /* perform FFT (in/out is data)*/
  fftw_execute_dft(fft.plan[1].our_fftw_plan, c_data, c_data);
  /* ===== second direction ===== */
  /* communication to current dir row format (in is data) */
  forw_grid_comm(fft.plan[2], data, fft.data_buf.data(), fft, comm);
  /* perform FFT (in/out is fft.data_buf)*/
  fftw_execute_dft(fft.plan[2].our_fftw_plan, c_data_buf, c_data_buf);
  /* ===== third direction  ===== */
  /* communication to current dir row format (in is fft.data_buf) */
  forw_grid_comm(fft.plan[3], fft.data_buf.data(), data, fft, comm);
  /* perform FFT (in/out is data)*/
  fftw_execute_dft(fft.plan[3].our_fftw_plan, c_data, c_data);

  /* REMARK: Result has to be in data. */
}

void fft_perform_back(double *data, bool check_complex, fft_data_struct &fft,
                      const boost::mpi::communicator &comm) {

  auto *c_data = (fftw_complex *)data;
  auto *c_data_buf = (fftw_complex *)fft.data_buf.data();

  /* ===== third direction  ===== */

  /* perform FFT (in is data) */
  fftw_execute_dft(fft.back[3].our_fftw_plan, c_data, c_data);
  /* communicate (in is data)*/
  back_grid_comm(fft.plan[3], fft.back[3], data, fft.data_buf.data(), fft,
                 comm);

  /* ===== second direction ===== */
  /* perform FFT (in is fft.data_buf) */
  fftw_execute_dft(fft.back[2].our_fftw_plan, c_data_buf, c_data_buf);
  /* communicate (in is fft.data_buf) */
  back_grid_comm(fft.plan[2], fft.back[2], fft.data_buf.data(), data, fft,
                 comm);

  /* ===== first direction  ===== */
  /* perform FFT (in is data) */
  fftw_execute_dft(fft.back[1].our_fftw_plan, c_data, c_data);
  /* throw away the (hopefully) empty complex component (in is data)*/
  for (int i = 0; i < fft.plan[1].new_size; i++) {
    fft.data_buf[i] = data[2 * i]; /* real value */
    // Vincent:
    if (check_complex && (data[2 * i + 1] > 1e-5)) {
      printf("Complex value is not zero (i=%d,data=%g)!!!\n", i,
             data[2 * i + 1]);
      if (i > 100)
        throw std::runtime_error("Complex value is not zero");
    }
  }
  /* communicate (in is fft.data_buf) */
  back_grid_comm(fft.plan[1], fft.back[1], fft.data_buf.data(), data, fft,
                 comm);

  /* REMARK: Result has to be in data. */
}

void fft_pack_block(double const *const in, double *const out,
                    int const start[3], int const size[3], int const dim[3],
                    int element) {
  /* linear index of in grid, linear index of out grid */
  int li_in, li_out = 0;
  /* copy size */
  int copy_size;
  /* offsets for indices in input grid */
  int m_in_offset, s_in_offset;
  /* offsets for indices in output grid */
  int m_out_offset;

  copy_size = element * size[2] * static_cast<int>(sizeof(double));
  m_in_offset = element * dim[2];
  s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  m_out_offset = element * size[2];
  li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (int s = 0; s < size[0]; s++) {
    for (int m = 0; m < size[1]; m++) {
      memmove(&(out[li_out]), &(in[li_in]), copy_size);
      li_in += m_in_offset;
      li_out += m_out_offset;
    }
    li_in += s_in_offset;
  }
}

void fft_unpack_block(double const *const in, double *const out,
                      int const start[3], int const size[3], int const dim[3],
                      int element) {
  /* mid and slow changing indices */
  int m, s;
  /* linear index of in grid, linear index of out grid */
  int li_in = 0, li_out;
  /* copy size */
  int copy_size;
  /* offset for indices in input grid */
  int m_in_offset;
  /* offsets for indices in output grid */
  int m_out_offset, s_out_offset;

  copy_size = element * size[2] * static_cast<int>(sizeof(double));
  m_out_offset = element * dim[2];
  s_out_offset = element * (dim[2] * (dim[1] - size[1]));
  m_in_offset = element * size[2];
  li_out = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (s = 0; s < size[0]; s++) {
    for (m = 0; m < size[1]; m++) {
      memmove(&(out[li_out]), &(in[li_in]), copy_size);
      li_in += m_in_offset;
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}
#endif
