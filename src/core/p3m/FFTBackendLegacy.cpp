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

#include "config/config.hpp"

#if defined(P3M) or defined(DP3M)

#include "FFTBackendLegacy.hpp"

#include "communication.hpp"

#include "fft/fft.hpp"

#include <memory>

template <typename FloatType>
FFTBackendLegacy<FloatType>::FFTBackendLegacy(P3MLocalMesh const &local_mesh)
    : FFTBackend<FloatType>(local_mesh),
      fft{std::make_unique<fft::fft_data_struct<FloatType>>(
          ::Communication::mpiCallbacksHandle()->share_mpi_env())} {}

template <typename FloatType>
FFTBackendLegacy<FloatType>::~FFTBackendLegacy() = default;

template <typename FloatType>
void FFTBackendLegacy<FloatType>::init(P3MParameters const &params) {
  ca_mesh_size = fft->initialize_fft(
      ::comm_cart, local_mesh.dim, local_mesh.margin, params.mesh,
      params.mesh_off, ks_pnum, ::communicator.node_grid);
}

template <typename FloatType>
void FFTBackendLegacy<FloatType>::forward_fft(FloatType *rs_mesh) {
  fft->forward_fft(::comm_cart, rs_mesh);
}

template <typename FloatType>
void FFTBackendLegacy<FloatType>::backward_fft(FloatType *rs_mesh) {
  fft->backward_fft(::comm_cart, rs_mesh, check_complex_residuals);
}

template class FFTBackendLegacy<float>;
template class FFTBackendLegacy<double>;

#endif // defined(P3M) or defined(DP3M)
