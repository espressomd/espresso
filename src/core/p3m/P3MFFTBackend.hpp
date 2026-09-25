/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#include <utils/Vector.hpp>

#include <complex>
#include <type_traits>

/**
 * @brief Abstract interface for the P3M reciprocal-space FFT.
 *
 * P3M only ever touches the FFT through this handful of methods: the k-space
 * box geometry accessors and the forward/backward transforms. Concrete
 * backends (heFFTe for the general/distributed case, kokkos-fft for the
 * single-rank case) implement it, so the solver is agnostic to which one is
 * active. The transform buffers are flat, contiguous, no-halo arrays; halo
 * exchange is handled outside the FFT.
 */
template <typename FloatType, class FFTConfig> struct P3MFFTBackend {
  using ComplexType = std::complex<FloatType>;
  /**
   * Real-space scalar type. P3M's real-space buffers (charge density, E-field)
   * are physically real in both configs, so this is always @c FloatType. For a
   * c2c backend (@c use_r2c=false) the imaginary part is dropped on the
   * backward transform via heFFTe's real-output convenience overload, matching
   * the pre-backend-interface behaviour where the templated forward/backward
   * transforms of @c P3MFFT accepted real pointers for either config.
   */
  using RSpaceScalar = FloatType;

  virtual ~P3MFFTBackend() = default;

  /** @brief Lower-left corner of this rank's k-space box (global coords). */
  virtual Utils::Vector3i ks_local_ld_index() const = 0;
  /** @brief Upper-right corner (exclusive) of this rank's k-space box. */
  virtual Utils::Vector3i ks_local_ur_index() const = 0;
  /** @brief Extent of this rank's k-space box. */
  virtual Utils::Vector3i ks_local_size() const = 0;
  /** @brief Extent of this rank's real-space (no-halo) box. */
  virtual Utils::Vector3i rs_local_size() const = 0;
  /** @brief Persistent, transform-ready buffer for the caller to fill with the
   *  no-halo real-space input; holds product(rs_local_size()) elements. Passing
   *  it back as @ref forward's @c in lets a backend transform it in place. */
  virtual RSpaceScalar *forward_input_buffer() = 0;
  /** @brief Forward transform (real/complex charge density to k-space). */
  virtual void forward(RSpaceScalar const *in, ComplexType *out) = 0;
  /** @brief Backward transform (k-space field to real space).
   *  The input buffer is non-const by contract: a backend may destroy it
   *  during the transform (FFTW's multi-dimensional c2r has no
   *  preserve-input mode, so the kokkos-fft backend runs it in place).
   *  Callers must recompute the k-space input before reading it again.
   */
  virtual void backward(ComplexType *in, RSpaceScalar *out) = 0;
};
