/*
 * Copyright (C) 2010-2024 The ESPResSo project
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
#include "vector.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>
#include <utils/math/permute_ifield.hpp>

#include <boost/mpi/communicator.hpp>

#include <fftw3.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

using Utils::get_linear_index;
using Utils::permute_ifield;

/** @name MPI tags for FFT communication */
/**@{*/
/** Tag for communication in forw_grid_comm() */
#define REQ_FFT_FORW 301
/** Tag for communication in back_grid_comm() */
#define REQ_FFT_BACK 302
/**@}*/

namespace fft {

/** This ugly function does the bookkeeping: which nodes have to
 *  communicate to each other, when you change the node grid.
 *  Changing the regular decomposition requires communication. This
 *  function finds (hopefully) the best way to do this. As input it
 *  needs the two grids (@p grid1, @p grid2) and a linear list (@p node_list1)
 *  with the node identities for @p grid1. The linear list (@p node_list2)
 *  for the second grid is calculated. For the communication group of
 *  the calling node it calculates a list (@c group) with the node
 *  identities and the positions (@p my_pos, @p pos) of that nodes in @p grid1
 *  and @p grid2. The return value is the size of the communication
 *  group. It gives -1 if the two grids do not fit to each other
 *  (@p grid1 and @p grid2 have to be component-wise multiples of each
 *  other, see e.g. \ref calc_2d_grid for how to do this).
 *
 *  \param[in]  grid1       The node grid you start with.
 *  \param[in]  grid2       The node grid you want to have.
 *  \param[in]  node_list1  Linear node index list for grid1.
 *  \param[out] node_list2  Linear node index list for grid2.
 *  \param[out] pos         Positions of the nodes in grid2
 *  \param[out] my_pos      Position of comm.rank() in grid2.
 *  \param[in]  rank        MPI rank.
 *  \return Size of the communication group.
 */
std::optional<std::vector<int>>
find_comm_groups(Utils::Vector3i const &grid1, Utils::Vector3i const &grid2,
                 std::span<int const> node_list1, std::span<int> node_list2,
                 std::span<int> pos, std::span<int> my_pos, int rank) {
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
  if (Utils::product(grid1) != Utils::product(grid2))
    return std::nullopt; /* unlike number of nodes */
  for (i = 0; i < 3; i++) {
    s1[i] = grid1[i] / grid2[i];
    if (s1[i] == 0)
      s1[i] = 1;
    else if (grid1[i] != grid2[i] * s1[i])
      return std::nullopt; /* grids do not match!!! */

    s2[i] = grid2[i] / grid1[i];
    if (s2[i] == 0)
      s2[i] = 1;
    else if (grid2[i] != grid1[i] * s2[i])
      return std::nullopt; /* grids do not match!!! */

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
          if (n == rank && my_group == 0) {
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
  return {group};
}

namespace {
/** Calculate the local fft mesh. Calculate the local mesh (@p loc_mesh)
 *  of a node at position (@p n_pos) in a node grid (@p n_grid) for a global
 *  mesh of size (@p mesh) and a mesh offset (@p mesh_off (in mesh units))
 *  and store also the first point (@p start) of the local mesh.
 *
 * \param[in]  n_pos    Position of the node in @p n_grid.
 * \param[in]  n_grid   node grid.
 * \param[in]  mesh     global mesh dimensions.
 * \param[in]  mesh_off global mesh offset.
 * \param[out] loc_mesh local mesh dimension.
 * \param[out] start    first point of local mesh in global mesh.
 * \return Number of mesh points in local mesh.
 */
int calc_local_mesh(const int *n_pos, const int *n_grid, const int *mesh,
                    const double *mesh_off, int *loc_mesh, int *start) {
  int last[3], size = 1;

  for (int i = 0; i < 3; i++) {
    auto const ai = mesh[i] / static_cast<double>(n_grid[i]);
    start[i] = static_cast<int>(ceil(ai * n_pos[i] - mesh_off[i]));
    last[i] = static_cast<int>(floor(ai * (n_pos[i] + 1) - mesh_off[i]));
    /* correct round off errors */
    if (ai * (n_pos[i] + 1) - mesh_off[i] - last[i] < 1.0e-15)
      last[i]--;
    if (1.0 + ai * n_pos[i] - mesh_off[i] - start[i] < 1.0e-15)
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
 *  \param[in]  mesh_off global mesh offset.
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

  /* offsets for indices in input grid */
  auto const m_in_offset = element * (dim[2] - size[2]);
  auto const s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  /* offset for mid changing indices of output grid */
  auto const m_out_offset = (element * size[0]) - element;
  /* linear index of in grid */
  int li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (int s = 0; s < size[0]; s++) { /* fast changing out */
    /* linear index of out grid */
    int li_out = element * s;
    for (int m = 0; m < size[1]; m++) {   /* slow changing out */
      for (int f = 0; f < size[2]; f++) { /* mid  changing out */
        for (int e = 0; e < element; e++)
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

  /* offsets for indices in input grid */
  auto const m_in_offset = element * (dim[2] - size[2]);
  auto const s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  /* offset for slow changing index of output grid */
  auto const s_out_offset = (element * size[0] * size[1]) - element;
  /* linear index of in grid */
  int li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (int s = 0; s < size[0]; s++) { /* mid changing out */
    auto const m_out_start = element * (s * size[1]);
    for (int m = 0; m < size[1]; m++) { /* fast changing out */
      /* linear index of out grid */
      int li_out = m_out_start + element * m;
      for (int f = 0; f < size[2]; f++) { /* slow  changing out */
        for (int e = 0; e < element; e++)
          out[li_out++] = in[li_in++];
        li_out += s_out_offset;
      }
      li_in += m_in_offset;
    }
    li_in += s_in_offset;
  }
}

} // namespace

/** Communicate the grid data according to the given forward FFT plan.
 *  \param comm   MPI communicator.
 *  \param plan   FFT communication plan.
 *  \param in     input mesh.
 *  \param out    output mesh.
 */
void fft_data_struct::forw_grid_comm(boost::mpi::communicator const &comm,
                                     fft_forw_plan const &plan,
                                     double const *in, double *out) {
  for (std::size_t i = 0ul; i < plan.group.size(); i++) {
    plan.pack_function(in, send_buf.data(), &(plan.send_block[6ul * i]),
                       &(plan.send_block[6ul * i + 3ul]), plan.old_mesh,
                       plan.element);

    if (plan.group[i] != comm.rank()) {
      MPI_Sendrecv(send_buf.data(), plan.send_size[i], MPI_DOUBLE,
                   plan.group[i], REQ_FFT_FORW, recv_buf.data(),
                   plan.recv_size[i], MPI_DOUBLE, plan.group[i], REQ_FFT_FORW,
                   comm, MPI_STATUS_IGNORE);
    } else { /* Self communication... */
      std::swap(send_buf, recv_buf);
    }
    fft_unpack_block(recv_buf.data(), out, &(plan.recv_block[6ul * i]),
                     &(plan.recv_block[6ul * i + 3ul]), plan.new_mesh,
                     plan.element);
  }
}

/** Communicate the grid data according to the given backward FFT plan.
 *  \param comm   MPI communicator.
 *  \param plan_f Forward FFT plan.
 *  \param plan_b Backward FFT plan.
 *  \param in     input mesh.
 *  \param out    output mesh.
 */
void fft_data_struct::back_grid_comm(boost::mpi::communicator const &comm,
                                     fft_forw_plan const &plan_f,
                                     fft_back_plan const &plan_b,
                                     double const *in, double *out) {
  /* Back means: Use the send/receive stuff from the forward plan but
     replace the receive blocks by the send blocks and vice
     versa. Attention then also new_mesh and old_mesh are exchanged */

  for (std::size_t i = 0ul; i < plan_f.group.size(); i++) {
    plan_b.pack_function(in, send_buf.data(), &(plan_f.recv_block[6ul * i]),
                         &(plan_f.recv_block[6ul * i + 3ul]), plan_f.new_mesh,
                         plan_f.element);

    if (plan_f.group[i] != comm.rank()) { /* send first, receive second */
      MPI_Sendrecv(send_buf.data(), plan_f.recv_size[i], MPI_DOUBLE,
                   plan_f.group[i], REQ_FFT_BACK, recv_buf.data(),
                   plan_f.send_size[i], MPI_DOUBLE, plan_f.group[i],
                   REQ_FFT_BACK, comm, MPI_STATUS_IGNORE);
    } else { /* Self communication... */
      std::swap(send_buf, recv_buf);
    }
    fft_unpack_block(recv_buf.data(), out, &(plan_f.send_block[6ul * i]),
                     &(plan_f.send_block[6ul * i + 3ul]), plan_f.old_mesh,
                     plan_f.element);
  }
}

/** Calculate 'best' mapping between a 2D and 3D grid.
 *  Required for the communication from 3D regular domain
 *  decomposition to 2D regular row decomposition.
 *  The dimensions of the 2D grid are resorted, if necessary, in a way
 *  that they are multiples of the 3D grid dimensions.
 *  \param g3d      3D grid.
 *  \param g2d      2D grid.
 *  \return         index of the row direction [0,1,2].
 */
int map_3don2d_grid(int const g3d[3], int g2d[3]) {
  int row_dir = -1;
  /* trivial case */
  if (g3d[2] == 1) {
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
  return row_dir;
}

/** Calculate most square 2D grid. */
static void calc_2d_grid(int n, int grid[3]) {
  grid[0] = n;
  grid[1] = 1;
  grid[2] = 1;
  for (auto i = static_cast<int>(std::sqrt(n)); i >= 1; i--) {
    if (n % i == 0) {
      grid[0] = n / i;
      grid[1] = i;
      grid[2] = 1;
      return;
    }
  }
}

int fft_data_struct::initialize_fft(boost::mpi::communicator const &comm,
                                    Utils::Vector3i const &ca_mesh_dim,
                                    int const *ca_mesh_margin,
                                    Utils::Vector3i const &global_mesh_dim,
                                    Utils::Vector3d const &global_mesh_off,
                                    int &ks_pnum, Utils::Vector3i const &grid) {

  int n_grid[4][3];         /* The four node grids. */
  int my_pos[4][3];         /* The position of comm.rank() in the node grids. */
  std::vector<int> n_id[4]; /* linear node identity lists for the node grids. */
  std::vector<int> n_pos[4]; /* positions of nodes in the node grids. */

  int const rank = comm.rank();
  int node_pos[3];
  MPI_Cart_coords(comm, rank, 3, node_pos);

  max_comm_size = 0;
  max_mesh_size = 0;
  for (int i = 0; i < 4; i++) {
    n_id[i].resize(1 * comm.size());
    n_pos[i].resize(3 * comm.size());
  }

  /* === node grids === */
  /* real space node grid (n_grid[0]) */
  for (int i = 0; i < 3; i++) {
    n_grid[0][i] = grid[i];
    my_pos[0][i] = node_pos[i];
  }
  for (int i = 0; i < comm.size(); i++) {
    MPI_Cart_coords(comm, i, 3, &(n_pos[0][3 * i + 0]));
    auto const lin_ind = get_linear_index(
        n_pos[0][3 * i + 0], n_pos[0][3 * i + 1], n_pos[0][3 * i + 2],
        {n_grid[0][0], n_grid[0][1], n_grid[0][2]});
    n_id[0][lin_ind] = i;
  }

  /* FFT node grids (n_grid[1 - 3]) */
  calc_2d_grid(comm.size(), n_grid[1]);
  /* resort n_grid[1] dimensions if necessary */
  forw[1].row_dir = map_3don2d_grid(n_grid[0], n_grid[1]);
  forw[0].n_permute = 0;
  for (int i = 1; i < 4; i++)
    forw[i].n_permute = (forw[1].row_dir + i) % 3;
  for (int i = 0; i < 3; i++) {
    n_grid[2][i] = n_grid[1][(i + 1) % 3];
    n_grid[3][i] = n_grid[1][(i + 2) % 3];
  }
  forw[2].row_dir = (forw[1].row_dir - 1) % 3;
  forw[3].row_dir = (forw[1].row_dir - 2) % 3;

  /* === communication groups === */
  /* copy local mesh off real space charge assignment grid */
  for (int i = 0; i < 3; i++)
    forw[0].new_mesh[i] = ca_mesh_dim[i];

  for (int i = 1; i < 4; i++) {
    auto group = find_comm_groups(
        {n_grid[i - 1][0], n_grid[i - 1][1], n_grid[i - 1][2]},
        {n_grid[i][0], n_grid[i][1], n_grid[i][2]}, n_id[i - 1],
        std::span(n_id[i]), std::span(n_pos[i]), my_pos[i], rank);
    if (not group) {
      /* try permutation */
      std::swap(n_grid[i][(forw[i].row_dir + 1) % 3],
                n_grid[i][(forw[i].row_dir + 2) % 3]);

      group = find_comm_groups(
          {n_grid[i - 1][0], n_grid[i - 1][1], n_grid[i - 1][2]},
          {n_grid[i][0], n_grid[i][1], n_grid[i][2]}, std::span(n_id[i - 1]),
          std::span(n_id[i]), std::span(n_pos[i]), my_pos[i], rank);

      if (not group) {
        throw std::runtime_error("INTERNAL ERROR: fft_find_comm_groups error");
      }
    }

    forw[i].group = group.value();

    forw[i].send_block.resize(6 * forw[i].group.size());
    forw[i].send_size.resize(forw[i].group.size());
    forw[i].recv_block.resize(6 * forw[i].group.size());
    forw[i].recv_size.resize(forw[i].group.size());

    forw[i].new_size = calc_local_mesh(
        my_pos[i], n_grid[i], global_mesh_dim.data(), global_mesh_off.data(),
        forw[i].new_mesh, forw[i].start);
    permute_ifield(forw[i].new_mesh, 3, -(forw[i].n_permute));
    permute_ifield(forw[i].start, 3, -(forw[i].n_permute));
    forw[i].n_ffts = forw[i].new_mesh[0] * forw[i].new_mesh[1];

    /* === send/recv block specifications === */
    for (std::size_t j = 0ul; j < forw[i].group.size(); j++) {
      /* send block: comm.rank() to comm-group-node i (identity: node) */
      int node = forw[i].group[j];
      forw[i].send_size[j] = calc_send_block(
          my_pos[i - 1], n_grid[i - 1], &(n_pos[i][3 * node]), n_grid[i],
          global_mesh_dim.data(), global_mesh_off.data(),
          &(forw[i].send_block[6ul * j]));
      permute_ifield(&(forw[i].send_block[6ul * j]), 3,
                     -(forw[i - 1].n_permute));
      permute_ifield(&(forw[i].send_block[6ul * j + 3ul]), 3,
                     -(forw[i - 1].n_permute));
      if (forw[i].send_size[j] > max_comm_size)
        max_comm_size = forw[i].send_size[j];
      /* First plan send blocks have to be adjusted, since the CA grid
         may have an additional margin outside the actual domain of the
         node */
      if (i == 1) {
        for (std::size_t k = 0ul; k < 3ul; k++)
          forw[1].send_block[6ul * j + k] += ca_mesh_margin[2ul * k];
      }
      /* recv block: comm.rank() from comm-group-node i (identity: node) */
      forw[i].recv_size[j] = calc_send_block(
          my_pos[i], n_grid[i], &(n_pos[i - 1][3 * node]), n_grid[i - 1],
          global_mesh_dim.data(), global_mesh_off.data(),
          &(forw[i].recv_block[6ul * j]));
      permute_ifield(&(forw[i].recv_block[6ul * j]), 3, -(forw[i].n_permute));
      permute_ifield(&(forw[i].recv_block[6ul * j + 3ul]), 3,
                     -(forw[i].n_permute));
      if (forw[i].recv_size[j] > max_comm_size)
        max_comm_size = forw[i].recv_size[j];
    }

    for (std::size_t j = 0ul; j < 3ul; j++)
      forw[i].old_mesh[j] = forw[i - 1].new_mesh[j];
    if (i == 1) {
      forw[i].element = 1;
    } else {
      forw[i].element = 2;
      for (std::size_t j = 0ul; j < forw[i].group.size(); j++) {
        forw[i].send_size[j] *= 2;
        forw[i].recv_size[j] *= 2;
      }
    }
  }

  /* Factor 2 for complex fields */
  max_comm_size *= 2;
  max_mesh_size = Utils::product(ca_mesh_dim);
  for (int i = 1; i < 4; i++)
    if (2 * forw[i].new_size > max_mesh_size)
      max_mesh_size = 2 * forw[i].new_size;

  /* === pack function === */
  for (int i = 1; i < 4; i++) {
    forw[i].pack_function = pack_block_permute2;
  }
  ks_pnum = 6;
  if (forw[1].row_dir == 2) {
    forw[1].pack_function = fft_pack_block;
    ks_pnum = 4;
  } else if (forw[1].row_dir == 1) {
    forw[1].pack_function = pack_block_permute1;
    ks_pnum = 5;
  }

  send_buf.resize(max_comm_size);
  recv_buf.resize(max_comm_size);
  data_buf.resize(max_mesh_size);
  auto *c_data = (fftw_complex *)(data_buf.data());

  /* === FFT Routines (Using FFTW / RFFTW package)=== */
  for (int i = 1; i < 4; i++) {
    if (init_tag) {
      forw[i].destroy_plan();
    }
    forw[i].dir = FFTW_FORWARD;
    forw[i].plan_handle =
        fftw_plan_many_dft(1, &forw[i].new_mesh[2], forw[i].n_ffts, c_data,
                           nullptr, 1, forw[i].new_mesh[2], c_data, nullptr, 1,
                           forw[i].new_mesh[2], forw[i].dir, FFTW_PATIENT);
    assert(forw[i].plan_handle);
  }

  /* === The BACK Direction === */
  /* this is needed because slightly different functions are used */
  for (int i = 1; i < 4; i++) {
    if (init_tag) {
      back[i].destroy_plan();
    }
    back[i].dir = FFTW_BACKWARD;
    back[i].plan_handle =
        fftw_plan_many_dft(1, &forw[i].new_mesh[2], forw[i].n_ffts, c_data,
                           nullptr, 1, forw[i].new_mesh[2], c_data, nullptr, 1,
                           forw[i].new_mesh[2], back[i].dir, FFTW_PATIENT);
    back[i].pack_function = pack_block_permute1;
    assert(back[i].plan_handle);
  }
  if (forw[1].row_dir == 2) {
    back[1].pack_function = fft_pack_block;
  } else if (forw[1].row_dir == 1) {
    back[1].pack_function = pack_block_permute2;
  }

  init_tag = true;

  return max_mesh_size;
}

void fft_data_struct::forward_fft(boost::mpi::communicator const &comm,
                                  double *data) {
  /* ===== first direction  ===== */

  auto *c_data = (fftw_complex *)data;
  auto *c_data_buf = (fftw_complex *)data_buf.data();

  /* communication to current dir row format (in is data) */
  forw_grid_comm(comm, forw[1], data, data_buf.data());

  /* complexify the real data array (in is data_buf) */
  for (int i = 0; i < forw[1].new_size; i++) {
    data[2 * i + 0] = data_buf[i]; /* real value */
    data[2 * i + 1] = 0;           /* complex value */
  }
  /* perform FFT (in/out is data)*/
  fftw_execute_dft(forw[1].plan_handle, c_data, c_data);
  /* ===== second direction ===== */
  /* communication to current dir row format (in is data) */
  forw_grid_comm(comm, forw[2], data, data_buf.data());
  /* perform FFT (in/out is data_buf) */
  fftw_execute_dft(forw[2].plan_handle, c_data_buf, c_data_buf);
  /* ===== third direction  ===== */
  /* communication to current dir row format (in is data_buf) */
  forw_grid_comm(comm, forw[3], data_buf.data(), data);
  /* perform FFT (in/out is data)*/
  fftw_execute_dft(forw[3].plan_handle, c_data, c_data);

  /* REMARK: Result has to be in data. */
}

void fft_data_struct::backward_fft(boost::mpi::communicator const &comm,
                                   double *data, bool check_complex) {

  auto *c_data = (fftw_complex *)data;
  auto *c_data_buf = (fftw_complex *)data_buf.data();

  /* ===== third direction  ===== */

  /* perform FFT (in is data) */
  fftw_execute_dft(back[3].plan_handle, c_data, c_data);
  /* communicate (in is data)*/
  back_grid_comm(comm, forw[3], back[3], data, data_buf.data());

  /* ===== second direction ===== */
  /* perform FFT (in is data_buf) */
  fftw_execute_dft(back[2].plan_handle, c_data_buf, c_data_buf);
  /* communicate (in is data_buf) */
  back_grid_comm(comm, forw[2], back[2], data_buf.data(), data);

  /* ===== first direction  ===== */
  /* perform FFT (in is data) */
  fftw_execute_dft(back[1].plan_handle, c_data, c_data);
  /* throw away the (hopefully) empty complex component (in is data) */
  for (int i = 0; i < forw[1].new_size; i++) {
    data_buf[i] = data[2 * i]; /* real value */
    // Vincent:
    if (check_complex and std::abs(data[2 * i + 1]) > 1e-5) {
      printf("Complex value is not zero (i=%d,data=%g)!!!\n", i,
             data[2 * i + 1]);
      if (i > 100)
        throw std::runtime_error("Complex value is not zero");
    }
  }
  /* communicate (in is data_buf) */
  back_grid_comm(comm, forw[1], back[1], data_buf.data(), data);

  /* REMARK: Result has to be in data. */
}

void fft_pack_block(double const *const in, double *const out,
                    int const start[3], int const size[3], int const dim[3],
                    int element) {

  auto const copy_size = element * size[2] * static_cast<int>(sizeof(double));
  /* offsets for indices in input grid */
  auto const m_in_offset = element * dim[2];
  auto const s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  /* offsets for indices in output grid */
  auto const m_out_offset = element * size[2];
  /* linear index of in grid, linear index of out grid */
  int li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));
  int li_out = 0;

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

  auto const copy_size = element * size[2] * static_cast<int>(sizeof(double));
  /* offsets for indices in output grid */
  auto const m_out_offset = element * dim[2];
  auto const s_out_offset = element * (dim[2] * (dim[1] - size[1]));
  /* offset for indices in input grid */
  auto const m_in_offset = element * size[2];
  /* linear index of in grid, linear index of out grid */
  int li_in = 0;
  int li_out = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (int s = 0; s < size[0]; s++) {
    for (int m = 0; m < size[1]; m++) {
      memmove(&(out[li_out]), &(in[li_in]), copy_size);
      li_in += m_in_offset;
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}

void fft_plan::destroy_plan() {
  if (plan_handle) {
    fftw_destroy_plan(plan_handle);
    plan_handle = nullptr;
  }
}

namespace detail {
void fft_free(void *const p) { ::fftw_free(p); }
void *fft_malloc(std::size_t length) { return ::fftw_malloc(length); }
} // namespace detail
} // namespace fft
