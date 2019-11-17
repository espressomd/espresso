//
// Created by florian on 17.11.19.
//

#include "p3m_send_mesh.hpp"
#include "fft.hpp"

#include <utils/mpi/cart_comm.hpp>

void p3m_send_mesh::resize(const boost::mpi::communicator &comm,
                           const p3m_local_mesh &local_mesh) {
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
  max = 0;
  for (int i = 0; i < 6; i++) {
    s_size[i] = 1;
    for (int j = 0; j < 3; j++) {
      s_dim[i][j] = s_ur[i][j] - s_ld[i][j];
      s_size[i] *= s_dim[i][j];
    }
    if (s_size[i] > max)
      max = s_size[i];
  }
  /* communication */
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(comm);

  int r_margin[6];
  for (int i = 0; i < 6; i++) {
    auto const j = (i % 2 == 0) ? i + 1 : i - 1;

    if (node_neighbors[i] != comm.rank()) {
      MPI_Sendrecv(&(local_mesh.margin[i]), 1, MPI_INT, node_neighbors[i],
                   REQ_P3M_INIT, &(r_margin[j]), 1, MPI_INT, node_neighbors[j],
                   REQ_P3M_INIT, comm, MPI_STATUS_IGNORE);
    } else {
      r_margin[j] = local_mesh.margin[i];
    }
  }
  /* recv grids */
  for (int i = 0; i < 3; i++)
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
  for (int i = 0; i < 6; i++) {
    r_size[i] = 1;
    for (int j = 0; j < 3; j++) {
      r_dim[i][j] = r_ur[i][j] - r_ld[i][j];
      r_size[i] *= r_dim[i][j];
    }
    if (r_size[i] > max)
      max = r_size[i];
  }

  send_grid.resize(max);
  recv_grid.resize(max);
}

void p3m_send_mesh::gather_grid(double *themesh,
                                const boost::mpi::communicator &comm,
                                const int dim[3]) {
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(comm);

  /* direction loop */
  for (int s_dir = 0; s_dir < 6; s_dir++) {
    auto const r_dir = (s_dir % 2 == 0) ? s_dir + 1 : s_dir - 1;

    /* pack send block */
    if (s_size[s_dir] > 0)
      fft_pack_block(themesh, send_grid.data(), s_ld[s_dir], s_dim[s_dir], dim,
                     1);

    /* communication */
    if (node_neighbors[s_dir] != comm.rank()) {
      MPI_Sendrecv(send_grid.data(), s_size[s_dir], MPI_DOUBLE,
                   node_neighbors[s_dir], REQ_P3M_GATHER, recv_grid.data(),
                   r_size[r_dir], MPI_DOUBLE, node_neighbors[r_dir],
                   REQ_P3M_GATHER, comm, MPI_STATUS_IGNORE);
    } else {
      std::swap(send_grid, recv_grid);
    }
    /* add recv block */
    if (r_size[r_dir] > 0) {
      p3m_add_block(recv_grid.data(), themesh, r_ld[r_dir], r_dim[r_dir], dim);
    }
  }
}

void p3m_send_mesh::spread_grid(double *themesh,
                                const boost::mpi::communicator &comm,
                                const int dim[3]) {
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(comm);

  /* direction loop */
  for (int s_dir = 5; s_dir >= 0; s_dir--) {
    auto const r_dir = (s_dir % 2 == 0) ? s_dir + 1 : s_dir - 1;

    /* pack send block */
    if (s_size[s_dir] > 0)
      fft_pack_block(themesh, send_grid.data(), r_ld[r_dir], r_dim[r_dir], dim,
                     1);
    /* communication */
    if (node_neighbors[r_dir] != comm.rank()) {
      MPI_Sendrecv(send_grid.data(), r_size[r_dir], MPI_DOUBLE,
                   node_neighbors[r_dir], REQ_P3M_SPREAD, recv_grid.data(),
                   s_size[s_dir], MPI_DOUBLE, node_neighbors[s_dir],
                   REQ_P3M_SPREAD, comm, MPI_STATUS_IGNORE);
    } else {
      std::swap(send_grid, recv_grid);
    }
    /* un pack recv block */
    if (s_size[s_dir] > 0) {
      fft_unpack_block(recv_grid.data(), themesh, s_ld[s_dir], s_dim[s_dir],
                       dim, 1);
    }
  }
}
