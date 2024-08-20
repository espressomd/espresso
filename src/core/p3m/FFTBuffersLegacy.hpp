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
#include <type_traits>

/** @brief Buffers for @ref FFTBackendLegacy. */
template <typename FloatType>
class FFTBuffersLegacy : public FFTBuffers<FloatType> {
  static_assert(std::is_same_v<FloatType, float> or
                    std::is_same_v<FloatType, double>,
                "FFTW only implements float and double");
  using FFTBuffers<FloatType>::FFTBuffers;
  using FFTBuffers<FloatType>::local_mesh;
  p3m_send_mesh<FloatType> mesh_comm;
  /** @brief real-space mesh (local) for CA/FFT. */
  fft::vector<FloatType> rs_mesh;
  /** @brief real-space mesh (local) for the electric or dipolar field. */
  std::array<fft::vector<FloatType>, 3u> rs_mesh_fields;

public:
  ~FFTBuffersLegacy() override;
  void init_halo() override;
  void init_meshes(int ca_mesh_size) override;
  void perform_scalar_halo_gather() override;
  void perform_vector_halo_gather() override;
  void perform_scalar_halo_spread() override;
  void perform_vector_halo_spread() override;
  void update_mesh_views(P3MFFTMesh<FloatType> &out) override;
  FloatType *get_scalar_mesh() override;
  std::array<FloatType *, 3u> get_vector_mesh() override;
};

#endif // defined(P3M) or defined(DP3M)
