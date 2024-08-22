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

#include "FFTBuffersLegacy.hpp"

#include "communication.hpp"

#include <array>
#include <span>

template <typename FloatType>
FFTBuffersLegacy<FloatType>::~FFTBuffersLegacy() = default;

template <typename FloatType>
void FFTBuffersLegacy<FloatType>::update_mesh_views(
    P3MFFTMesh<FloatType> &out) {
  out.rs_scalar = std::span(rs_mesh);
  for (auto i = 0u; i < 3u; ++i) {
    out.rs_fields[i] = std::span(rs_mesh_fields[i]);
  }
}

template <typename FloatType> void FFTBuffersLegacy<FloatType>::init_halo() {
  mesh_comm.resize(::comm_cart, local_mesh);
}

template <typename FloatType>
void FFTBuffersLegacy<FloatType>::init_meshes(int ca_mesh_size) {
  rs_mesh.resize(ca_mesh_size);
  for (auto &rs_mesh_field : rs_mesh_fields) {
    rs_mesh_field.resize(ca_mesh_size);
  }
}

template <typename FloatType>
void FFTBuffersLegacy<FloatType>::perform_vector_halo_spread() {
  std::array<FloatType *, 3u> meshes = {{rs_mesh_fields[0u].data(),
                                         rs_mesh_fields[1u].data(),
                                         rs_mesh_fields[2u].data()}};
  mesh_comm.spread_grid(::comm_cart, meshes, local_mesh.dim);
}

template <typename FloatType>
void FFTBuffersLegacy<FloatType>::perform_scalar_halo_gather() {
  mesh_comm.gather_grid(::comm_cart, rs_mesh.data(), local_mesh.dim);
}

template <typename FloatType>
void FFTBuffersLegacy<FloatType>::perform_vector_halo_gather() {
  std::array<FloatType *, 3u> meshes = {{rs_mesh_fields[0u].data(),
                                         rs_mesh_fields[1u].data(),
                                         rs_mesh_fields[2u].data()}};
  mesh_comm.gather_grid(::comm_cart, meshes, local_mesh.dim);
}

template <typename FloatType>
void FFTBuffersLegacy<FloatType>::perform_scalar_halo_spread() {
  mesh_comm.spread_grid(::comm_cart, rs_mesh.data(), local_mesh.dim);
}

template <typename FloatType>
FloatType *FFTBuffersLegacy<FloatType>::get_scalar_mesh() {
  return rs_mesh.data();
}

template <typename FloatType>
std::array<FloatType *, 3u> FFTBuffersLegacy<FloatType>::get_vector_mesh() {
  return {rs_mesh_fields[0u].data(), rs_mesh_fields[1u].data(),
          rs_mesh_fields[2u].data()};
}

template class FFTBuffersLegacy<float>;
template class FFTBuffersLegacy<double>;

#endif // defined(P3M) or defined(DP3M)
