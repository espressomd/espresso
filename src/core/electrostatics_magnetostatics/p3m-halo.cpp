#include "p3m-halo.hpp"
#include "pack_block.hpp"

#include <utils/mpi/cart_comm.hpp>

#include <boost/range/numeric.hpp>

p3m_halo_comm calc_send_mesh(const boost::mpi::communicator &comm,
                             const int dim[3], const int margin[6]) {
  p3m_halo_comm send_mesh;

  send_mesh.comm = comm;

  int done[3] = {0, 0, 0};
  /* send grids */
  for (int i = 0; i < 3; i++) {
    send_mesh.dim[i] = dim[i];

    for (int j = 0; j < 3; j++) {
      /* left */
      send_mesh.s_ld[i * 2][j] = 0 + done[j] * margin[j * 2];
      if (j == i)
        send_mesh.s_ur[i * 2][j] = margin[j * 2];
      else
        send_mesh.s_ur[i * 2][j] = dim[j] - done[j] * margin[(j * 2) + 1];
      /* right */
      if (j == i)
        send_mesh.s_ld[(i * 2) + 1][j] = (dim[j] - margin[j * 2 + 1]);
      else
        send_mesh.s_ld[(i * 2) + 1][j] = 0 + done[j] * margin[j * 2];
      send_mesh.s_ur[(i * 2) + 1][j] = dim[j] - done[j] * margin[(j * 2) + 1];
    }
    done[i] = 1;
  }
  send_mesh.max = 0;
  for (int i = 0; i < 6; i++) {
    send_mesh.s_size[i] = 1;
    for (int j = 0; j < 3; j++) {
      send_mesh.s_dim[i][j] = send_mesh.s_ur[i][j] - send_mesh.s_ld[i][j];
      send_mesh.s_size[i] *= send_mesh.s_dim[i][j];
    }
    if (send_mesh.s_size[i] > send_mesh.max)
      send_mesh.max = send_mesh.s_size[i];
  }
  /* communication */
  auto const node_neighbors = Utils::Mpi::calc_face_neighbors<3>(comm);

  int r_margin[6];
  for (int i = 0; i < 6; i++) {
    auto const j = (i % 2 == 0) ? i + 1 : i - 1;

    MPI_Sendrecv(&(margin[i]), 1, MPI_INT, node_neighbors[i], 200,
                 &(r_margin[j]), 1, MPI_INT, node_neighbors[j], 200, comm,
                 MPI_STATUS_IGNORE);
  }
  /* recv grids */
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      if (j == i) {
        send_mesh.r_ld[i * 2][j] = send_mesh.s_ld[i * 2][j] + margin[2 * j];
        send_mesh.r_ur[i * 2][j] = send_mesh.s_ur[i * 2][j] + r_margin[2 * j];
        send_mesh.r_ld[(i * 2) + 1][j] =
            send_mesh.s_ld[(i * 2) + 1][j] - r_margin[(2 * j) + 1];
        send_mesh.r_ur[(i * 2) + 1][j] =
            send_mesh.s_ur[(i * 2) + 1][j] - margin[(2 * j) + 1];
      } else {
        send_mesh.r_ld[i * 2][j] = send_mesh.s_ld[i * 2][j];
        send_mesh.r_ur[i * 2][j] = send_mesh.s_ur[i * 2][j];
        send_mesh.r_ld[(i * 2) + 1][j] = send_mesh.s_ld[(i * 2) + 1][j];
        send_mesh.r_ur[(i * 2) + 1][j] = send_mesh.s_ur[(i * 2) + 1][j];
      }
    }
  for (int i = 0; i < 6; i++) {
    send_mesh.r_size[i] = 1;
    for (int j = 0; j < 3; j++) {
      send_mesh.r_dim[i][j] = send_mesh.r_ur[i][j] - send_mesh.r_ld[i][j];
      send_mesh.r_size[i] *= send_mesh.r_dim[i][j];
    }
    if (send_mesh.r_size[i] > send_mesh.max)
      send_mesh.max = send_mesh.r_size[i];
  }

  return send_mesh;
}

void p3m_gather_halo(Utils::Span<double *const> data,
                     const p3m_halo_comm &send_mesh) {
  auto const node_neighbors =
      Utils::Mpi::calc_face_neighbors<3>(send_mesh.comm);

  auto const buf_size = send_mesh.max * data.size();
  send_mesh.send_buffer.resize(buf_size);
  send_mesh.recv_buffer.resize(buf_size);

  /* direction loop */
  for (int s_dir = 0; s_dir < 6; s_dir++) {
    auto const r_dir = (s_dir % 2 == 0) ? s_dir + 1 : s_dir - 1;

    /* pack send block */
    if (send_mesh.s_size[s_dir] > 0) {
      boost::accumulate(data, send_mesh.send_buffer.data(),
                        [&](double *send_buf, double const *in_buf) {
                          return pack_block(
                              in_buf, send_buf, send_mesh.s_ld[s_dir],
                              send_mesh.s_dim[s_dir], send_mesh.dim, 1,
                              Utils::MemoryOrder::ROW_MAJOR);
                        });
    }

    /* communication */
    if (node_neighbors[s_dir] != send_mesh.comm.rank()) {
      MPI_Sendrecv(
          send_mesh.send_buffer.data(), data.size() * send_mesh.s_size[s_dir],
          MPI_DOUBLE, node_neighbors[s_dir], 201, send_mesh.recv_buffer.data(),
          data.size() * send_mesh.r_size[r_dir], MPI_DOUBLE,
          node_neighbors[r_dir], 201, send_mesh.comm, MPI_STATUS_IGNORE);
    } else {
      std::swap(send_mesh.send_buffer, send_mesh.recv_buffer);
    }

    /* add recv blocks */
    if (send_mesh.r_size[r_dir] > 0) {
      boost::accumulate(
          data, static_cast<const double *>(send_mesh.recv_buffer.data()),
          [&](const double *recv_buf, double *out_buf) {
            return unpack_block(recv_buf, out_buf, send_mesh.r_ld[r_dir],
                                send_mesh.r_dim[r_dir], send_mesh.dim, 1,
                                Utils::MemoryOrder::ROW_MAJOR, std::plus<>());
          });
    }
  }
}

void p3m_gather_halo(double *data, const p3m_halo_comm &send_mesh) {
  p3m_gather_halo(Utils::make_const_span(&data, 1), send_mesh);
}

void p3m_spread_halo(Utils::Span<double *const> data,
                     const p3m_halo_comm &send_mesh) {
  auto const node_neighbors =
      Utils::Mpi::calc_face_neighbors<3>(send_mesh.comm);

  /* Make sure the buffers are large enough */
  auto const buf_size = send_mesh.max * data.size();
  send_mesh.send_buffer.resize(buf_size);
  send_mesh.recv_buffer.resize(buf_size);

  /* direction loop */
  for (int s_dir = 5; s_dir >= 0; s_dir--) {
    auto const r_dir = (s_dir % 2 == 0) ? s_dir + 1 : s_dir - 1;

    /* pack send block */
    if (send_mesh.s_size[s_dir] > 0) {
      boost::accumulate(data, send_mesh.send_buffer.data(),
                        [&](double *send_buf, const double *in_buf) {
                          return pack_block(
                              in_buf, send_buf, send_mesh.r_ld[r_dir],
                              send_mesh.r_dim[r_dir], send_mesh.dim, 1,
                              Utils::MemoryOrder::ROW_MAJOR);
                        });
    }

    /* communication */
    if (node_neighbors[r_dir] != send_mesh.comm.rank()) {
      MPI_Sendrecv(
          send_mesh.send_buffer.data(), data.size() * send_mesh.r_size[r_dir],
          MPI_DOUBLE, node_neighbors[r_dir], 202, send_mesh.recv_buffer.data(),
          data.size() * send_mesh.s_size[s_dir], MPI_DOUBLE,
          node_neighbors[s_dir], 202, send_mesh.comm, MPI_STATUS_IGNORE);
    } else {
      std::swap(send_mesh.recv_buffer, send_mesh.send_buffer);
    }
    /* un pack recv block */
    if (send_mesh.s_size[s_dir] > 0) {
      boost::accumulate(
          data, static_cast<const double *>(send_mesh.recv_buffer.data()),
          [&](const double *recv_buf, double *out_buf) {
            return unpack_block(recv_buf, out_buf, send_mesh.s_ld[s_dir],
                                send_mesh.s_dim[s_dir], send_mesh.dim, 1,
                                Utils::MemoryOrder::ROW_MAJOR);
          });
    }
  }
}

void p3m_spread_halo(double *data, const p3m_halo_comm &send_mesh) {
  p3m_spread_halo(Utils::make_const_span(&data, 1), send_mesh);
}
