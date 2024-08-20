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

#pragma once

/** \file
 *
 *  Routines, row decomposition, data structures and communication for the
 *  3D-FFT.
 *
 *  The 3D-FFT is split into three 1D-FFTs. The data is
 *  distributed in such a way, that for the actual direction of the
 *  FFT each node has a certain number of rows for which it performs a
 *  1D-FFT. After performing the FFT on that direction the data is
 *  redistributed.
 *
 *  For simplicity, a full complex-to-complex FFT is implemented,
 *  even though a real-to-complex FFT would be sufficient.
 */

#include "vector.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

struct fftw_plan_s;
struct fftwf_plan_s;
namespace boost::mpi {
class environment;
class communicator;
} // namespace boost::mpi

namespace fft {

template <typename FloatType> struct fft_plan {
  static_assert(std::is_same_v<FloatType, float> or
                    std::is_same_v<FloatType, double>,
                "FFTW only implements float and double");
  using fftw_plan = std::conditional_t<std::is_same_v<FloatType, double>,
                                       fftw_plan_s *, fftwf_plan_s *>;

  ~fft_plan() { destroy_plan(); }

  /** plan direction: forward or backward FFT (enum value from FFTW). */
  int dir;
  /** plan for the FFT. */
  fftw_plan plan_handle = nullptr;
  /** packing function for send blocks. */
  void (*pack_function)(FloatType const *const, FloatType *const, int const *,
                        int const *, int const *, int);
  void destroy_plan();
};

/** @brief Plan for a forward 1D FFT of a flattened 3D array. */
template <typename FloatType>
struct fft_forw_plan : public fft_plan<FloatType> {
  /** row direction of that FFT. */
  int row_dir;
  /** permutations from normal coordinate system. */
  int n_permute;
  /** number of 1D FFTs. */
  int n_ffts;

  /** size of local mesh before communication. */
  std::array<int, 3u> old_mesh;
  /** size of local mesh after communication, also used for actual FFT. */
  std::array<int, 3u> new_mesh;
  /** lower left point of local FFT mesh in global FFT mesh coordinates. */
  std::array<int, 3u> start;
  /** size of new mesh (number of mesh points). */
  int new_size;

  /** group of nodes which have to communicate with each other. */
  std::vector<int> group;

  /** Send block specification. 6 integers for each node: start[3], size[3]. */
  std::vector<int> send_block;
  /** Send block communication sizes. */
  std::vector<int> send_size;
  /** Recv block specification. 6 integers for each node: start[3], size[3]. */
  std::vector<int> recv_block;
  /** Recv block communication sizes. */
  std::vector<int> recv_size;
  /** size of send block elements, i.e. 1 for real, 2 for complex. */
  int element;
};

/** @brief Plan for a backward 1D FFT of a flattened 3D array. */
template <typename FloatType>
struct fft_back_plan : public fft_plan<FloatType> {};

/**
 * @brief Information about the three one dimensional FFTs and how the nodes
 * have to communicate inbetween.
 *
 * @note FFT numbering starts with 1 for technical reasons (because we have 4
 * node grids, the index 0 is used for the real space charge assignment grid).
 */
template <typename FloatType> struct fft_data_struct {
private:
  /**
   * @brief Handle to the MPI environment.
   * Has to be the first member in the class definition, so that FFT plans
   * are destroyed before the MPI environment expires (non-static class
   * members are destroyed in the reverse order of their initialization).
   */
  std::shared_ptr<boost::mpi::environment> m_mpi_env;

  /** Information for forward FFTs. */
  std::array<fft_forw_plan<FloatType>, 4u> forw;
  /** Information for backward FFTs. */
  std::array<fft_back_plan<FloatType>, 4u> back;

  /** Whether FFT is initialized or not. */
  bool init_tag = false;

  /** Maximal size of the communication buffers. */
  int max_comm_size = 0;

  /** Maximal local mesh size. */
  int max_mesh_size = 0;

  /** send buffer. */
  std::vector<FloatType> send_buf;
  /** receive buffer. */
  std::vector<FloatType> recv_buf;
  /** Buffer for receive data. */
  fft::vector<FloatType> data_buf;

public:
  explicit fft_data_struct(decltype(m_mpi_env) mpi_env)
      : m_mpi_env{std::move(mpi_env)} {}

  // disable copy construction: unsafe because we store raw pointers
  // to FFT plans (avoids double-free and use-after-free)
  fft_data_struct &operator=(fft_data_struct<FloatType> const &) = delete;
  fft_data_struct(fft_data_struct<FloatType> const &) = delete;

  /** Initialize everything connected to the 3D-FFT.
   *
   *  \param[in]  comm            MPI communicator.
   *  \param[in]  ca_mesh_dim     Local CA mesh dimensions.
   *  \param[in]  ca_mesh_margin  Local CA mesh margins.
   *  \param[in]  global_mesh_dim Global CA mesh dimensions.
   *  \param[in]  global_mesh_off Global CA mesh offset.
   *  \param[out] ks_pnum         Number of permutations in k-space.
   *  \param[in]  grid            Number of nodes in each spatial dimension.
   *  \return Maximal size of local fft mesh (needed for allocation of ca_mesh).
   */
  int initialize_fft(boost::mpi::communicator const &comm,
                     Utils::Vector3i const &ca_mesh_dim,
                     int const *ca_mesh_margin,
                     Utils::Vector3i const &global_mesh_dim,
                     Utils::Vector3d const &global_mesh_off, int &ks_pnum,
                     Utils::Vector3i const &grid);

  /** Perform an in-place forward 3D FFT.
   *  \warning The content of \a data is overwritten.
   *  \param[in,out] data  Mesh.
   *  \param[in]     comm  MPI communicator
   */
  void forward_fft(boost::mpi::communicator const &comm, FloatType *data);

  /** Perform an in-place backward 3D FFT.
   *  \warning The content of \a data is overwritten.
   *  \param[in,out] data           Mesh.
   *  \param[in]     check_complex  Throw an error if the complex component is
   *                                non-zero.
   *  \param[in]     comm           MPI communicator.
   */
  void backward_fft(boost::mpi::communicator const &comm, FloatType *data,
                    bool check_complex);

  auto const &get_mesh_size() const { return forw[3u].new_mesh; }

  auto const &get_mesh_start() const { return forw[3u].start; }

private:
  void forw_grid_comm(boost::mpi::communicator const &comm,
                      fft_forw_plan<FloatType> const &plan, FloatType const *in,
                      FloatType *out);
  void back_grid_comm(boost::mpi::communicator const &comm,
                      fft_forw_plan<FloatType> const &plan_f,
                      fft_back_plan<FloatType> const &plan_b,
                      FloatType const *in, FloatType *out);
};

int map_3don2d_grid(int const g3d[3], int g2d[3]);

std::optional<std::vector<int>> find_comm_groups(Utils::Vector3i const &,
                                                 Utils::Vector3i const &,
                                                 std::span<int const>,
                                                 std::span<int>, std::span<int>,
                                                 std::span<int>, int);

} // namespace fft
