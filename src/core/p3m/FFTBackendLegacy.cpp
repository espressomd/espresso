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

#include <utils/Vector.hpp>

#include <array>
#include <memory>
#include <span>
#include <utility>

template <typename FloatType>
FFTBackendLegacy<FloatType>::FFTBackendLegacy(
    p3m_data_struct_fft<FloatType> &obj, bool dipolar)
    : FFTBackend<FloatType>(obj), dipolar{dipolar},
      fft{std::make_unique<fft::fft_data_struct<FloatType>>(
          ::Communication::mpiCallbacksHandle()->share_mpi_env())} {}

template <typename FloatType>
FFTBackendLegacy<FloatType>::~FFTBackendLegacy() = default;

template <typename FloatType>
void FFTBackendLegacy<FloatType>::update_mesh_data() {
  auto const mesh_size_ptr = fft->get_mesh_size();
  auto const mesh_start_ptr = fft->get_mesh_start();
  for (auto i = 0u; i < 3u; ++i) {
    mesh.size[i] = mesh_size_ptr[i];
    mesh.start[i] = mesh_start_ptr[i];
  }
  mesh.stop = mesh.start + mesh.size;
  mesh.ks_scalar = std::span(ks_mesh);
  mesh.rs_scalar = std::span(rs_mesh);
  for (auto i = 0u; i < 3u; ++i) {
    mesh.rs_fields[i] = std::span(rs_mesh_fields[i]);
  }
}

template <typename FloatType> void FFTBackendLegacy<FloatType>::init_fft() {
  mesh_comm.resize(::comm_cart, local_mesh);
  auto ca_mesh_size = fft->initialize_fft(
      ::comm_cart, local_mesh.dim, local_mesh.margin, params.mesh,
      params.mesh_off, mesh.ks_pnum, ::communicator.node_grid);
  rs_mesh.resize(ca_mesh_size);
  if (dipolar) {
    ks_mesh.resize(ca_mesh_size);
  }
  for (auto &rs_mesh_field : rs_mesh_fields) {
    rs_mesh_field.resize(ca_mesh_size);
  }
  update_mesh_data();
}

template <typename FloatType>
void FFTBackendLegacy<FloatType>::perform_vector_back_fft() {
  /* Back FFT force component mesh */
  for (auto &rs_mesh_field : rs_mesh_fields) {
    fft->backward_fft(::comm_cart, rs_mesh_field.data(),
                      check_complex_residuals);
  }
  /* redistribute force component mesh */
  std::array<FloatType *, 3u> meshes = {{rs_mesh_fields[0u].data(),
                                         rs_mesh_fields[1u].data(),
                                         rs_mesh_fields[2u].data()}};
  mesh_comm.spread_grid(::comm_cart, meshes, local_mesh.dim);
}

template <typename FloatType>
void FFTBackendLegacy<FloatType>::perform_scalar_fwd_fft() {
  mesh_comm.gather_grid(::comm_cart, rs_mesh.data(), local_mesh.dim);
  fft->forward_fft(::comm_cart, rs_mesh.data());
  update_mesh_data();
}

template <typename FloatType>
void FFTBackendLegacy<FloatType>::perform_vector_fwd_fft() {
  std::array<FloatType *, 3u> meshes = {{rs_mesh_fields[0u].data(),
                                         rs_mesh_fields[1u].data(),
                                         rs_mesh_fields[2u].data()}};
  mesh_comm.gather_grid(::comm_cart, meshes, local_mesh.dim);
  for (auto &rs_mesh_field : rs_mesh_fields) {
    fft->forward_fft(::comm_cart, rs_mesh_field.data());
  }
  update_mesh_data();
}

template <typename FloatType>
void FFTBackendLegacy<FloatType>::perform_scalar_back_fft() {
  /* Back FFT force component mesh */
  fft->backward_fft(::comm_cart, rs_mesh.data(), check_complex_residuals);
  /* redistribute force component mesh */
  mesh_comm.spread_grid(::comm_cart, rs_mesh.data(), local_mesh.dim);
}

template class FFTBackendLegacy<float>;
template class FFTBackendLegacy<double>;

#endif // defined(P3M) or defined(DP3M)
