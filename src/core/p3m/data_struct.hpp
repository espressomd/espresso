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

#pragma once

#include <config/config.hpp>

#if defined(ESPRESSO_P3M) or defined(ESPRESSO_DP3M)

#include "common.hpp"
#include "communication.hpp"

#include <utils/device_qualifier.hpp>

#include <array>
#include <cassert>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

/** @brief State of the p3m methods, the part which applies to  both,
  electrostatic and dipolar p3m */
template <typename FloatType, Arch Architecture = Arch::CPU>
struct P3MStateCommon {
  using value_type = FloatType;

  explicit P3MStateCommon(P3MParameters &&parameters)
      : params{std::move(parameters)}, local_mesh{} {}

  /** @brief P3M base parameters. */
  P3MParameters params;
  /** @brief Local mesh geometry information for this MPI rank */
  P3MLocalMesh local_mesh;

  /**
   * @brief Spatial differential operator in k-space.
   * We use an i*k differentiation.
   * TODO: RW: I think this is not the differential operator but the mapping
   * between index in the GLOBAL mesh and the corresponding
   * k-vector with a few pre-factors missing.
   */
  std::array<std::vector<int>, 3> d_op;

  /** Calculate the Fourier transformed differential operator.
   *  Remark: This is done on the level of n-vectors and not k-vectors,
   *  i.e. the prefactor @f$ 2i\pi/L @f$ is missing!
   */
  void calc_differential_operator() {
    d_op = calc_p3m_mesh_shift(params.mesh, true);
  }

  /** @brief Force optimised influence function (k-space) */
  std::vector<FloatType> g_force;
  /** @brief Energy optimised influence function (k-space) */
  std::vector<FloatType> g_energy;
};

#if defined(ESPRESSO_DP3M)

/**
 * @brief API for the legacy FFT backend of the P3M algorithm.
 */
template <typename FloatType> class FFTBackend {
protected:
  P3MLocalMesh const &local_mesh;

public:
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
  virtual std::array<int, 3u> const &get_mesh_size() const = 0;
  virtual std::array<int, 3u> const &get_mesh_start() const = 0;
};

/**
 * @brief API for the legacy FFT mesh buffers.
 */
template <typename FloatType> class FFTBuffers {
protected:
  P3MLocalMesh const &local_mesh;

public:
  explicit FFTBuffers(P3MLocalMesh const &local_mesh)
      : local_mesh{local_mesh} {}
  virtual HOST_ONLY_QUALIFIER ~FFTBuffers() = default;
  /** @brief Initialize the meshes. */
  virtual void init_meshes(int ca_mesh_size) = 0;
  /** @brief Initialize the halo buffers. */
  virtual void init_halo() = 0;
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

#endif // defined(ESPRESSO_DP3M)

#endif // defined(ESPRESSO_P3M) or defined(ESPRESSO_DP3M)
