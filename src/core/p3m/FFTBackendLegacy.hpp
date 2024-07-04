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
#include "send_mesh.hpp"

#include "fft/vector.hpp"

#include <array>
#include <memory>
#include <tuple>

namespace fft {
struct fft_data_struct;
} // namespace fft

/**
 * @brief Historic FFT backend based on FFTW3.
 * The 3D FFT is split into three 1D FFTs.
 */
class FFTBackendLegacy : public FFTBackend {
  bool dipolar;
  std::unique_ptr<fft::fft_data_struct> fft;
  /** @brief k-space mesh (local) for k-space calculations. */
  std::vector<double> ks_mesh;
  /** @brief real-space mesh (local) for CA/FFT. */
  fft::vector<double> rs_mesh;
  /** @brief real-space mesh (local) for the electric or dipolar field. */
  std::array<fft::vector<double>, 3> rs_mesh_fields;
  p3m_send_mesh mesh_comm;

public:
  explicit FFTBackendLegacy(p3m_data_struct &obj, bool dipolar);
  ~FFTBackendLegacy() override;
  void init_fft() override;
  void perform_fwd_fft() override;
  void perform_field_back_fft() override;
  void perform_space_back_fft() override;
  void update_mesh_data();

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
