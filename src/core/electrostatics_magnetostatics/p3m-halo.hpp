#ifndef ESPRESSO_P3M_HALO_HPP
#define ESPRESSO_P3M_HALO_HPP

#include <utils/Span.hpp>

#include <boost/mpi/communicator.hpp>

#include <utils/index.hpp>
#include <vector>

/** Structure for send/recv meshes. */
struct halo_comm {
  /** dimensions of the mesh */
  int dim[3];
  /** dimension of sub meshes to send. */
  int s_dim[6][3];
  /** left down corners of sub meshes to send. */
  int s_ld[6][3];
  /** up right corners of sub meshes to send. */
  int s_ur[6][3];
  /** sizes for send buffers. */
  int s_size[6];
  /** dimension of sub meshes to recv. */
  int r_dim[6][3];
  /** left down corners of sub meshes to recv. */
  int r_ld[6][3];
  /** up right corners of sub meshes to recv. */
  int r_ur[6][3];
  /** sizes for recv buffers. */
  int r_size[6];
  /** maximal size for send/recv buffers. */
  int max;

  /** MPI communicator for this halo comm */
  boost::mpi::communicator comm;

  /** vector to store grid points to send. */
  mutable std::vector<double> send_buffer;
  /** vector to store grid points to recv */
  mutable std::vector<double> recv_buffer;
};

halo_comm plan_halo_comm(const boost::mpi::communicator &comm,
                         const Utils::Vector3i &dim,
                         const Utils::Array<int, 6> &margin);

/**
 * @brief Add halo regions to their original images.
 * @param data The mesh data
 * @param send_mesh Halo plan
 */
void p3m_gather_halo(double *data, const halo_comm &send_mesh,
                     Utils::MemoryOrder memory_order);
void p3m_gather_halo(Utils::Span<double *const> data,
                     const halo_comm &send_mesh,
                     Utils::MemoryOrder memory_order);

/**
 * @brief Overwrite halo regions with their original images.
 * @param data The mesh data
 * @param send_mesh Halo plan
 */
void p3m_spread_halo(double *data, const halo_comm &send_mesh,
                     Utils::MemoryOrder memory_order);
void p3m_spread_halo(Utils::Span<double *const> data,
                     const halo_comm &send_mesh,
                     Utils::MemoryOrder memory_order);

#endif // ESPRESSO_P3M_HALO_HPP
