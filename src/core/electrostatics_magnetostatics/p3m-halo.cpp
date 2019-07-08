#include "p3m-halo.hpp"
#include "pack_block.hpp"

#include <utils/mpi/cart_comm.hpp>

#include <boost/range/numeric.hpp>

halo_comm::halo_comm(const boost::mpi::communicator &comm_,
                     const Utils::Vector3i &dim_,
                     const Utils::Array<int, 6> &margin) {
  comm = comm_;

  int done[3] = {0, 0, 0};
  /* send grids */
  for (int i = 0; i < 3; i++) {
    dim[i] = dim_[i];

    for (int j = 0; j < 3; j++) {
      /* left */
      s_ld[i * 2][j] = 0 + done[j] * margin[j * 2];
      if (j == i)
        s_ur[i * 2][j] = margin[j * 2];
      else
        s_ur[i * 2][j] = dim[j] - done[j] * margin[(j * 2) + 1];
      /* right */
      if (j == i)
        s_ld[(i * 2) + 1][j] = (dim[j] - margin[j * 2 + 1]);
      else
        s_ld[(i * 2) + 1][j] = 0 + done[j] * margin[j * 2];
      s_ur[(i * 2) + 1][j] = dim[j] - done[j] * margin[(j * 2) + 1];
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
        r_ld[i * 2][j] = s_ld[i * 2][j] + margin[2 * j];
        r_ur[i * 2][j] = s_ur[i * 2][j] + r_margin[2 * j];
        r_ld[(i * 2) + 1][j] = s_ld[(i * 2) + 1][j] - r_margin[(2 * j) + 1];
        r_ur[(i * 2) + 1][j] = s_ur[(i * 2) + 1][j] - margin[(2 * j) + 1];
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
}

void halo_comm::gather(Utils::Span<double *const> data,
                       Utils::MemoryOrder memory_order) const {
  auto const node_neighbors = Utils::Mpi::calc_face_neighbors<3>(comm);

  auto const buf_size = max * data.size();
  send_buffer.resize(buf_size);
  recv_buffer.resize(buf_size);

  /* direction loop */
  for (int s_dir = 0; s_dir < 6; s_dir++) {
    auto const r_dir = (s_dir % 2 == 0) ? s_dir + 1 : s_dir - 1;

    /* pack send block */
    if (s_size[s_dir] > 0) {
      boost::accumulate(data, send_buffer.data(),
                        [&](double *send_buf, double const *in_buf) {
                          return pack_block(in_buf, send_buf, s_ld[s_dir],
                                            s_dim[s_dir], dim, 1, memory_order);
                        });
    }

    /* communication */
    if (node_neighbors[s_dir] != comm.rank()) {
      MPI_Sendrecv(send_buffer.data(), data.size() * s_size[s_dir], MPI_DOUBLE,
                   node_neighbors[s_dir], 201, recv_buffer.data(),
                   data.size() * r_size[r_dir], MPI_DOUBLE,
                   node_neighbors[r_dir], 201, comm, MPI_STATUS_IGNORE);
    } else {
      std::swap(send_buffer, recv_buffer);
    }

    /* add recv blocks */
    if (r_size[r_dir] > 0) {
      boost::accumulate(data, static_cast<const double *>(recv_buffer.data()),
                        [&](const double *recv_buf, double *out_buf) {
                          return unpack_block(recv_buf, out_buf, r_ld[r_dir],
                                              r_dim[r_dir], dim, 1,
                                              memory_order, std::plus<>());
                        });
    }
  }
}

void halo_comm::gather(double *data, Utils::MemoryOrder memory_order) const {
  gather(Utils::make_const_span(&data, 1), memory_order);
}

void halo_comm::spread(Utils::Span<double *const> data,
                       Utils::MemoryOrder memory_order) const {
  auto const node_neighbors = Utils::Mpi::calc_face_neighbors<3>(comm);

  /* Make sure the buffers are large enough */
  auto const buf_size = max * data.size();
  send_buffer.resize(buf_size);
  recv_buffer.resize(buf_size);

  /* direction loop */
  for (int s_dir = 5; s_dir >= 0; s_dir--) {
    auto const r_dir = (s_dir % 2 == 0) ? s_dir + 1 : s_dir - 1;

    /* pack send block */
    if (s_size[s_dir] > 0) {
      boost::accumulate(data, send_buffer.data(),
                        [&](double *send_buf, const double *in_buf) {
                          return pack_block(in_buf, send_buf, r_ld[r_dir],
                                            r_dim[r_dir], dim, 1, memory_order);
                        });
    }

    /* communication */
    if (node_neighbors[r_dir] != comm.rank()) {
      MPI_Sendrecv(send_buffer.data(), data.size() * r_size[r_dir], MPI_DOUBLE,
                   node_neighbors[r_dir], 202, recv_buffer.data(),
                   data.size() * s_size[s_dir], MPI_DOUBLE,
                   node_neighbors[s_dir], 202, comm, MPI_STATUS_IGNORE);
    } else {
      std::swap(recv_buffer, send_buffer);
    }
    /* un pack recv block */
    if (s_size[s_dir] > 0) {
      boost::accumulate(data, static_cast<const double *>(recv_buffer.data()),
                        [&](const double *recv_buf, double *out_buf) {
                          return unpack_block(recv_buf, out_buf, s_ld[s_dir],
                                              s_dim[s_dir], dim, 1,
                                              memory_order);
                        });
    }
  }
}
