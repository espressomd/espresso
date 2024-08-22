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
#include "data_struct.hpp"

#include <array>
#include <memory>
#include <tuple>
#include <type_traits>

namespace fft {
template <typename FloatType> struct fft_data_struct;
} // namespace fft

/**
 * @brief Historic FFT backend based on FFTW3.
 * The 3D FFT is split into three 1D FFTs.
 */
template <typename FloatType>
class FFTBackendLegacy : public FFTBackend<FloatType> {
  static_assert(std::is_same_v<FloatType, float> or
                    std::is_same_v<FloatType, double>,
                "FFTW only implements float and double");
  std::unique_ptr<fft::fft_data_struct<FloatType>> fft;
  using FFTBackend<FloatType>::local_mesh;
  using FFTBackend<FloatType>::check_complex_residuals;
  int ca_mesh_size = -1;
  int ks_pnum = -1;

public:
  FFTBackendLegacy(P3MLocalMesh const &local_mesh);
  ~FFTBackendLegacy() override;
  void init(P3MParameters const &params) override;
  void forward_fft(FloatType *rs_mesh) override;
  void backward_fft(FloatType *rs_mesh) override;
  int get_ca_mesh_size() const noexcept override { return ca_mesh_size; }
  int get_ks_pnum() const noexcept override { return ks_pnum; }
  std::array<int, 3u> const &get_mesh_size() const override {
    return fft->get_mesh_size();
  }
  std::array<int, 3u> const &get_mesh_start() const override {
    return fft->get_mesh_start();
  }

  /**
   * @brief Index helpers for reciprocal space.
   * After the FFT the data is in order YZX, which
   * means that Y is the slowest changing index.
   */
  std::tuple<int, int, int> get_permutations() const override {
    constexpr static int KX = 2;
    constexpr static int KY = 0;
    constexpr static int KZ = 1;
    return {KX, KY, KZ};
  }
};

#endif // defined(P3M) or defined(DP3M)
