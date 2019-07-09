#ifndef ESPRESSO_P3M_HALO_HPP
#define ESPRESSO_P3M_HALO_HPP

#include <utils/Span.hpp>
#include <utils/index.hpp>

#include <boost/mpi/communicator.hpp>

#include <vector>

/**
 * @brief Halo communication plan for regular meshes.
 */
class HaloComm {
  /** dimensions of the mesh */
  int dim[3];
  /** dimension of sub meshes to send. */
  int s_dim[6][3];
  /** left down corners of sub meshes to send. */
  int s_ld[6][3];
  /** up right corners of sub meshes to send. */
  int s_ur[6][3];
  /** sizes for send buffers. */
  int s_size[6] = {};
  /** dimension of sub meshes to recv. */
  int r_dim[6][3];
  /** left down corners of sub meshes to recv. */
  int r_ld[6][3];
  /** up right corners of sub meshes to recv. */
  int r_ur[6][3];
  /** sizes for recv buffers. */
  int r_size[6] = {};
  /** maximal size for send/recv buffers. */
  int max = 0;

  /** MPI communicator for this halo comm */
  boost::mpi::communicator comm;

  /** vector to store grid points to send. */
  mutable std::vector<double> send_buffer;
  /** vector to store grid points to recv */
  mutable std::vector<double> recv_buffer;

public:
  HaloComm() = default;
  /**
   * @param comm Valid cartesian communicator
   * @param dim Local grid dimensions including halo.
   * @param margin Size of the halo in each of the face directions.
   */
  HaloComm(const boost::mpi::communicator &comm, const Utils::Vector3i &dim,
           const Utils::Array<int, 6> &margin);

  /**
   * @brief Overwrite halo regions with their original images.
   * @param data The mesh data
   * @param memory_order row- or column-major
   */
  void spread(double *data, Utils::MemoryOrder memory_order) const {
    spread(Utils::make_const_span(&data, 1), memory_order);
  }

  /**
   * @brief Overwrite halo regions on multiple meshes with their original
   * images.
   * @param data Multiple meshes to update
   * @param memory_order row- or column-major
   */
  void spread(Utils::Span<double *const> data,
              Utils::MemoryOrder memory_order) const;

  /**
   * @brief Add halo regions to their original images.
   * @param data The mesh data
   * @param memory_order row- or column-major
   */
  void gather(double *data, Utils::MemoryOrder memory_order) const;

  /**
   * @brief Add halo regions to their original images on multiple meshes.
   * @param data Multiple meshes to update
   * @param memory_order row- or column-major
   */
  void gather(Utils::Span<double *const> data,
              Utils::MemoryOrder memory_order) const;
};
#endif // ESPRESSO_P3M_HALO_HPP
