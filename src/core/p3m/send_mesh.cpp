/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#include <config/config.hpp>

#if defined(ESPRESSO_P3M) or defined(ESPRESSO_DP3M)

#include "fft/fft.hpp"
#include "p3m/common.hpp"
#include "p3m/packing.hpp"
#include "p3m/send_mesh.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/cart_comm.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/datatype.hpp>

#include <mpi.h>

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <initializer_list>
#include <span>
#include <utility>

template <typename T>
static void mesh_sendrecv(T const *const sendbuf, int scount, int dest,
                          T *const recvbuf, int rcount, int source,
                          boost::mpi::communicator const &comm, int tag) {
  auto const type = boost::mpi::get_mpi_datatype<T>(*sendbuf);
  MPI_Sendrecv(reinterpret_cast<void const *>(sendbuf), scount, type, dest, tag,
               reinterpret_cast<void *>(recvbuf), rcount, type, source, tag,
               comm, MPI_STATUS_IGNORE);
}

/** Add values of a 3d-grid input block (size[3]) to values of 3d-grid
 *  output array with dimension dim[3] at start position start[3].
 *
 *  \param in          Pointer to first element of input block data.
 *  \param out         Pointer to first element of output grid.
 *  \param start       Start position of block in output grid.
 *  \param size        Dimensions of the block
 *  \param dim         Dimensions of the output grid.
 */
template <typename FloatType>
static void p3m_add_block(FloatType const *in, FloatType *out,
                          int const start[3], int const size[3],
                          int const dim[3]) {
  /* fast, mid and slow changing indices */
  int f, m, s;
  /* linear index of in grid, linear index of out grid */
  int li_in = 0, li_out = 0;
  /* offsets for indices in output grid */
  int m_out_offset, s_out_offset;

  li_out = start[2] + (dim[2] * (start[1] + (dim[1] * start[0])));
  m_out_offset = dim[2] - size[2];
  s_out_offset = (dim[2] * (dim[1] - size[1]));

  for (s = 0; s < size[0]; s++) {
    for (m = 0; m < size[1]; m++) {
      for (f = 0; f < size[2]; f++) {
        out[li_out++] += in[li_in++];
      }
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}

template <typename FloatType>
void p3m_send_mesh<FloatType>::resize(boost::mpi::communicator const &comm,
                                      P3MLocalMesh const &local_mesh) {
  int done[3] = {0, 0, 0};
  /* send grids */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      /* left */
      s_ld[i * 2][j] = 0 + done[j] * local_mesh.margin[j * 2];
      if (j == i)
        s_ur[i * 2][j] = local_mesh.margin[j * 2];
      else
        s_ur[i * 2][j] =
            local_mesh.dim[j] - done[j] * local_mesh.margin[(j * 2) + 1];
      /* right */
      if (j == i)
        s_ld[(i * 2) + 1][j] = local_mesh.in_ur[j];
      else
        s_ld[(i * 2) + 1][j] = 0 + done[j] * local_mesh.margin[j * 2];
      s_ur[(i * 2) + 1][j] =
          local_mesh.dim[j] - done[j] * local_mesh.margin[(j * 2) + 1];
    }
    done[i] = 1;
  }
  for (int i = 0; i < 6; i++) {
    s_size[i] = 1;
    for (int j = 0; j < 3; j++) {
      s_dim[i][j] = s_ur[i][j] - s_ld[i][j];
      s_size[i] *= s_dim[i][j];
    }
  }
  /* communication */
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(comm);

  int r_margin[6];
  for (int i = 0; i < 6; i++) {
    auto const j = (i % 2 == 0) ? i + 1 : i - 1;

    if (node_neighbors[i] != comm.rank()) {
      mesh_sendrecv(&(local_mesh.margin[i]), 1, node_neighbors[i],
                    &(r_margin[j]), 1, node_neighbors[j], comm, REQ_P3M_INIT);
    } else {
      r_margin[j] = local_mesh.margin[i];
    }
  }
  /* recv grids */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (j == i) {
        r_ld[i * 2][j] = s_ld[i * 2][j] + local_mesh.margin[2 * j];
        r_ur[i * 2][j] = s_ur[i * 2][j] + r_margin[2 * j];
        r_ld[(i * 2) + 1][j] = s_ld[(i * 2) + 1][j] - r_margin[(2 * j) + 1];
        r_ur[(i * 2) + 1][j] =
            s_ur[(i * 2) + 1][j] - local_mesh.margin[(2 * j) + 1];
      } else {
        r_ld[i * 2][j] = s_ld[i * 2][j];
        r_ur[i * 2][j] = s_ur[i * 2][j];
        r_ld[(i * 2) + 1][j] = s_ld[(i * 2) + 1][j];
        r_ur[(i * 2) + 1][j] = s_ur[(i * 2) + 1][j];
      }
    }
  }
  for (int i = 0; i < 6; i++) {
    r_size[i] = 1;
    for (int j = 0; j < 3; j++) {
      r_dim[i][j] = r_ur[i][j] - r_ld[i][j];
      r_size[i] *= r_dim[i][j];
    }
  }
  /* both directions of a dimension share the send/recv buffers */
  max = 0;
  for (int d = 0; d < 3; ++d) {
    auto const l = 2 * d;
    auto const r = l + 1;
    max = std::max(max, s_size[l] + s_size[r]);
    max = std::max(max, r_size[l] + r_size[r]);
  }
}

template <typename FloatType>
void p3m_send_mesh<FloatType>::exchange(boost::mpi::communicator const &comm,
                                        int left, int right, int send_left,
                                        int send_right, int recv_right,
                                        int recv_left, int tag) {
  if (left == comm.rank()) {
    // 1 MPI rank along this direction: exchange with self
    std::swap(send_grid, recv_grid);
  } else if (left == right) {
    // 2 MPI ranks along this direction: send both blocks in a single message
    mesh_sendrecv(send_grid.data(), send_left + send_right, left,
                  recv_grid.data(), recv_right + recv_left, left, comm, tag);
  } else {
    // 3 or more MPI ranks along this direction: overlapping communication
    auto const type = boost::mpi::get_mpi_datatype(FloatType{});
    std::array<MPI_Request, 4> reqs;
    MPI_Irecv(recv_grid.data(), recv_right, type, right, tag, comm, &reqs[0]);
    MPI_Irecv(recv_grid.data() + recv_right, recv_left, type, left, tag, comm,
              &reqs[1]);
    MPI_Isend(send_grid.data(), send_left, type, left, tag, comm, &reqs[2]);
    MPI_Isend(send_grid.data() + send_left, send_right, type, right, tag, comm,
              &reqs[3]);
    MPI_Waitall(4, reqs.data(), MPI_STATUSES_IGNORE);
  }
}

template <typename FloatType>
void p3m_send_mesh<FloatType>::gather_grid(boost::mpi::communicator const &comm,
                                           std::span<FloatType *> meshes,
                                           Utils::Vector3i const &dim) {
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(comm);
  auto const n_meshes = static_cast<int>(meshes.size());
  send_grid.resize(static_cast<std::size_t>(max) * meshes.size());
  recv_grid.resize(static_cast<std::size_t>(max) * meshes.size());

  /* dimension loop: both directions of a dimension are independent, since
   * margins are packed and received blocks are added to inner layers */
  for (int d = 0; d < 3; d++) {
    auto const l = 2 * d;
    auto const r = l + 1;

    /* pack send blocks */
    auto *out = send_grid.data();
    for (auto const dir : {l, r}) {
      for (auto *mesh : meshes) {
        fft_pack_block(mesh, out, s_ld[dir], s_dim[dir], dim.data(), 1);
        out += s_size[dir];
      }
    }

    /* communication */
    exchange(comm, node_neighbors[l], node_neighbors[r], n_meshes * s_size[l],
             n_meshes * s_size[r], n_meshes * r_size[r], n_meshes * r_size[l],
             REQ_P3M_GATHER);

    /* add recv blocks */
    auto const *in = recv_grid.data();
    for (auto const dir : {r, l}) {
      for (auto *mesh : meshes) {
        p3m_add_block(in, mesh, r_ld[dir], r_dim[dir], dim.data());
        in += r_size[dir];
      }
    }
  }
}

template <typename FloatType>
void p3m_send_mesh<FloatType>::spread_grid(boost::mpi::communicator const &comm,
                                           std::span<FloatType *> meshes,
                                           Utils::Vector3i const &dim) {
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(comm);
  auto const n_meshes = static_cast<int>(meshes.size());
  send_grid.resize(static_cast<std::size_t>(max) * meshes.size());
  recv_grid.resize(static_cast<std::size_t>(max) * meshes.size());

  /* dimension loop: both directions of a dimension are independent, since
   * inner layers are packed and received blocks overwrite the margins */
  for (int d = 2; d >= 0; --d) {
    auto const l = 2 * d;
    auto const r = l + 1;

    /* pack send blocks */
    auto *out = send_grid.data();
    for (auto const dir : {l, r}) {
      for (auto *mesh : meshes) {
        fft_pack_block(mesh, out, r_ld[dir], r_dim[dir], dim.data(), 1);
        out += r_size[dir];
      }
    }

    /* communication */
    exchange(comm, node_neighbors[l], node_neighbors[r], n_meshes * r_size[l],
             n_meshes * r_size[r], n_meshes * s_size[r], n_meshes * s_size[l],
             REQ_P3M_SPREAD);

    /* unpack recv blocks */
    auto const *in = recv_grid.data();
    for (auto const dir : {r, l}) {
      for (auto *mesh : meshes) {
        fft_unpack_block(in, mesh, s_ld[dir], s_dim[dir], dim.data(), 1);
        in += s_size[dir];
      }
    }
  }
}

template class p3m_send_mesh<float>;
template class p3m_send_mesh<double>;

#endif // defined(ESPRESSO_P3M) or defined(ESPRESSO_DP3M)
