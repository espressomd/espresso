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

/**
 * @brief Base class for the electrostatics and magnetostatics P3M algorithms.
 * Contains a handle to the FFT backend, information about the local mesh,
 * the differential operator, and various buffers.
 */
struct p3m_data_struct {
  explicit p3m_data_struct(P3MParameters &&parameters)
      : params{std::move(parameters)} {}

  /** @brief P3M base parameters. */
  P3MParameters params;
  /** @brief Local mesh properties. */
  P3MLocalMesh local_mesh;

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
};

template <typename FloatType>
struct p3m_data_struct_fft : public p3m_data_struct {
  using p3m_data_struct::p3m_data_struct;
  using value_type = FloatType;
  /** @brief Local mesh FFT buffers. */
  P3MFFTMesh<FloatType> mesh;

  /** Force optimised influence function (k-space) */
  std::vector<FloatType> g_force;
  /** Energy optimised influence function (k-space) */
  std::vector<FloatType> g_energy;
  /** FFT backend. */
  std::unique_ptr<FFTBackend<FloatType>> fft;

  template <typename T, class... Args> void make_fft_instance(Args... args) {
    assert(fft == nullptr);
    fft = std::make_unique<T>(*this, args...);
  }
};

/**
 * @brief API for the FFT backend of the P3M algorithm.
 * Any FFT backend must implement this interface.
 * The backend can read some members of @ref p3m_data_struct
 * but can only modify the FFT buffers in @ref P3MFFTMesh.
 */
template <typename FloatType> class FFTBackend {
protected:
  P3MParameters const &params;
  P3MLocalMesh const &local_mesh;
  P3MFFTMesh<FloatType> &mesh;

public:
  bool check_complex_residuals = false;
  explicit FFTBackend(p3m_data_struct_fft<FloatType> &obj)
      : params{obj.params}, local_mesh{obj.local_mesh}, mesh{obj.mesh} {}
  virtual ~FFTBackend() = default;
  /** @brief Initialize the FFT plans and buffers. */
  virtual void init_fft() = 0;
  /** @brief Carry out the forward FFT of the scalar mesh. */
  virtual void perform_scalar_fwd_fft() = 0;
  /** @brief Carry out the forward FFT of the vector meshes. */
  virtual void perform_vector_fwd_fft() = 0;
  /** @brief Carry out the backward FFT of the scalar mesh. */
  virtual void perform_scalar_back_fft() = 0;
  /** @brief Carry out the backward FFT of the vector meshes. */
  virtual void perform_vector_back_fft() = 0;
  /** @brief Get indices of the k-space data layout. */
  virtual std::tuple<int, int, int> get_permutations() const = 0;
};

#endif // defined(P3M) or defined(DP3M)
