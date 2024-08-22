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

#include "config/config.hpp"

#if defined(P3M) or defined(DP3M)

#include "common.hpp"

#include <array>
#include <cassert>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

template <typename FloatType> class FFTBackend;
template <typename FloatType> class FFTBuffers;

/**
 * @brief Base class for the electrostatics and magnetostatics P3M algorithms.
 * Contains a handle to the FFT backend, information about the local mesh,
 * the differential operator, and various buffers.
 */
template <typename FloatType> struct p3m_data_struct {
  using value_type = FloatType;

  explicit p3m_data_struct(P3MParameters &&parameters)
      : params{std::move(parameters)} {}

  /** @brief P3M base parameters. */
  P3MParameters params;
  /** @brief Local mesh properties. */
  P3MLocalMesh local_mesh;
  /** @brief Local mesh FFT buffers. */
  P3MFFTMesh<FloatType> mesh;

  /**
   * @brief Spatial differential operator in k-space.
   * We use an i*k differentiation.
   */
  std::array<std::vector<int>, 3> d_op;

  /** Calculate the Fourier transformed differential operator.
   *  Remark: This is done on the level of n-vectors and not k-vectors,
   *  i.e. the prefactor @f$ 2i\pi/L @f$ is missing!
   */
  void calc_differential_operator() {
    d_op = detail::calc_meshift(params.mesh, true);
  }

  /** @brief Force optimised influence function (k-space) */
  std::vector<FloatType> g_force;
  /** @brief Energy optimised influence function (k-space) */
  std::vector<FloatType> g_energy;
  /** @brief FFT algorithm. */
  std::unique_ptr<FFTBackend<FloatType>> fft;
  /** @brief FFT buffers. */
  std::unique_ptr<FFTBuffers<FloatType>> fft_buffers;

  void init();

  void update_mesh_views() {
    auto const mesh_size_ptr = fft->get_mesh_size();
    auto const mesh_start_ptr = fft->get_mesh_start();
    for (auto i = 0u; i < 3u; ++i) {
      mesh.size[i] = mesh_size_ptr[i];
      mesh.start[i] = mesh_start_ptr[i];
    }
    mesh.stop = mesh.start + mesh.size;
    fft_buffers->update_mesh_views(mesh);
  }

  template <typename T, class... Args> void make_fft_instance(Args... args) {
    assert(fft == nullptr);
    fft = std::make_unique<T>(std::as_const(local_mesh), args...);
  }

  template <typename T, class... Args> void make_mesh_instance(Args... args) {
    assert(fft_buffers == nullptr);
    fft_buffers = std::make_unique<T>(std::as_const(local_mesh), args...);
  }
};

/**
 * @brief API for the FFT backend of the P3M algorithm.
 */
template <typename FloatType> class FFTBackend {
protected:
  P3MLocalMesh const &local_mesh;

public:
  bool check_complex_residuals = false;
  explicit FFTBackend(P3MLocalMesh const &local_mesh)
      : local_mesh{local_mesh} {}
  virtual ~FFTBackend() = default;
  virtual void init(P3MParameters const &params) = 0;
  virtual int get_ca_mesh_size() const noexcept = 0;
  virtual int get_ks_pnum() const noexcept = 0;
  /** @brief Carry out the forward FFT of the scalar mesh. */
  virtual void forward_fft(FloatType *rs_mesh) = 0;
  /** @brief Carry out the backward FFT of the scalar mesh. */
  virtual void backward_fft(FloatType *rs_mesh) = 0;
  /** @brief Get indices of the k-space data layout. */
  virtual std::tuple<int, int, int> get_permutations() const = 0;
  virtual std::array<int, 3u> const &get_mesh_size() const = 0;
  virtual std::array<int, 3u> const &get_mesh_start() const = 0;
};

/**
 * @brief API for the FFT mesh buffers.
 */
template <typename FloatType> class FFTBuffers {
protected:
  P3MLocalMesh const &local_mesh;

public:
  bool check_complex_residuals = false;
  explicit FFTBuffers(P3MLocalMesh const &local_mesh)
      : local_mesh{local_mesh} {}
  virtual ~FFTBuffers() = default;
  /** @brief Initialize the meshes. */
  virtual void init_meshes(int ca_mesh_size) = 0;
  /** @brief Initialize the halo buffers. */
  virtual void init_halo() = 0;
  /** @brief Update scalar mesh halo with data from neighbors (accumulation). */
  virtual void perform_scalar_halo_gather() = 0;
  /** @brief Update vector mesh halo with data from neighbors (accumulation). */
  virtual void perform_vector_halo_gather() = 0;
  /** @brief Update scalar mesh halo of all neighbors. */
  virtual void perform_scalar_halo_spread() = 0;
  /** @brief Update vector mesh halo of all neighbors. */
  virtual void perform_vector_halo_spread() = 0;
  /**
   * @brief Get pointer to scalar mesh begin.
   * Should only be used by @ref FFTBackend.
   */
  virtual FloatType *get_scalar_mesh() = 0;
  /**
   * @brief Get pointer to vector mesh begin.
   * Should only be used by @ref FFTBackend.
   */
  virtual std::array<FloatType *, 3u> get_vector_mesh() = 0;
  /**
   * @brief Update the scalar and vector mesh views in @ref P3MFFTMesh
   * to point to the new underlying data structures.
   */
  virtual void update_mesh_views(P3MFFTMesh<FloatType> &out) = 0;
};

#endif // defined(P3M) or defined(DP3M)
