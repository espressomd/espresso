/*
 * Copyright (C) 2023-2026 The ESPResSo project
 * Copyright (C) 2020 The waLBerla project
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

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 007e77e077ad9d22b5eed6f3d3118240993e553c

/*
 * Lattice field accessors.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/a16141524c58ab88386e2a0f8fdd7c63c5edd704/python/lbmpy_walberla/templates/LatticeModel.tmpl.h
 */

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Vector3.h>

#include <field/iterators/IteratorMacros.h>

#include <gpu/FieldAccessor.h>
#include <gpu/FieldIndexing.h>
#include <gpu/GPUField.h>
#include <gpu/Kernel.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

#include <array>
#include <vector>

#if defined(__NVCC__)
#define RESTRICT __restrict__
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
#define RESTRICT __restrict__
#else
// clang compiling CUDA code in host mode
#define RESTRICT __restrict__
#endif
#endif
#elif defined(__GNUC__) or defined(__GNUG__)
#define RESTRICT __restrict__
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

/** @brief Get linear index of flattened data with original layout @c fzyx. */
static __forceinline__ __device__ uint getLinearIndex(uint3 blockIdx, uint3 threadIdx, uint3 gridDim, uint3 blockDim, uint fOffset) {
  auto const x = threadIdx.x;
  auto const y = blockIdx.x;
  auto const z = blockIdx.y;
  auto const f = blockIdx.z;
  auto const ySize = gridDim.x;
  auto const zSize = gridDim.y;
  auto const fSize = fOffset;
  return f +
         z * fSize +
         y * fSize * zSize +
         x * fSize * zSize * ySize;
}

namespace walberla {
namespace ek {
namespace accessor {

namespace Scalar {
// LCOV_EXCL_START
__global__ void kernel_get(
    gpu::FieldAccessor<float> scalar_field,
    float *out) {
  auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 1u);
  scalar_field.set(blockIdx, threadIdx);
  out += offset;
  if (scalar_field.isValidPosition()) {
    out[0u] = scalar_field.get(0u);
  }
}

__global__ void kernel_set(
    gpu::FieldAccessor<float> scalar_field,
    float const *RESTRICT in) {
  auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 1u);
  scalar_field.set(blockIdx, threadIdx);
  in += offset;
  if (scalar_field.isValidPosition()) {
    scalar_field.get(0u) = in[0u];
  }
}

__global__ void kernel_broadcast(
    gpu::FieldAccessor<float> scalar_field,
    float const in) {
  scalar_field.set(blockIdx, threadIdx);
  if (scalar_field.isValidPosition()) {
    scalar_field.get(0u) = in;
  }
}

__global__ void kernel_add(
    gpu::FieldAccessor<float> scalar_field,
    float const *RESTRICT in) {
  auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 1u);
  scalar_field.set(blockIdx, threadIdx);
  in += offset;
  if (scalar_field.isValidPosition()) {
    scalar_field.get(0u) += in[0u];
  }
}

__global__ void kernel_broadcast_add(
    gpu::FieldAccessor<float> scalar_field,
    float const in) {
  scalar_field.set(blockIdx, threadIdx);
  if (scalar_field.isValidPosition()) {
    scalar_field.get(0u) += in;
  }
}
// LCOV_EXCL_STOP

float get(
    gpu::GPUField<float> const *scalar_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<float> dev_data(1u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*scalar_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  float result{};
  thrust::copy(dev_data.begin(), dev_data.end(), &result);
  return result;
}

void set(
    gpu::GPUField<float> *scalar_field,
    float const value,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  auto kernel = gpu::make_kernel(kernel_broadcast);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*scalar_field, ci));
  kernel.addParam(value);
  kernel();
}

void add(
    gpu::GPUField<float> *scalar_field,
    float const value,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  auto kernel = gpu::make_kernel(kernel_add);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*scalar_field, ci));
  kernel.addParam(value);
  kernel();
}

void initialize(
    gpu::GPUField<float> *scalar_field,
    float const value) {
  CellInterval ci = scalar_field->xyzSizeWithGhostLayer();
  auto kernel = gpu::make_kernel(kernel_broadcast);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*scalar_field, ci));
  kernel.addParam(value);
  kernel();
}

void add_to_all(
    gpu::GPUField<float> *scalar_field,
    Vector3<float> const value) {
  CellInterval ci = scalar_field->xyzSizeWithGhostLayer();
  auto kernel = gpu::make_kernel(kernel_broadcast_add);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*scalar_field, ci));
  kernel.addParam(value);
  kernel();
}

std::vector<float> get(
    gpu::GPUField<float> const *scalar_field,
    CellInterval const &ci) {
  thrust::device_vector<float> dev_data(ci.numCells() * 1u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*scalar_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<float> out(ci.numCells());
  thrust::copy(dev_data.begin(), dev_data.end(), out.data());
  return out;
}

void set(
    gpu::GPUField<float> *scalar_field,
    std::vector<float> const &values,
    CellInterval const &ci) {
  thrust::device_vector<float> dev_data(values.begin(), values.end());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*scalar_field, ci));
  kernel.addParam(const_cast<const float *>(dev_data_ptr));
  kernel();
}
} // namespace Scalar

namespace Vector {
// LCOV_EXCL_START
__global__ void kernel_get(
    gpu::FieldAccessor<float> vec,
    float *u_out) {
  auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 3u);
  vec.set(blockIdx, threadIdx);
  u_out += offset;
  if (vec.isValidPosition()) {
    u_out[0u] = vec.get(0u);
    u_out[1u] = vec.get(1u);
    u_out[2u] = vec.get(2u);
  }
}

__global__ void kernel_set(
    gpu::FieldAccessor<float> vec,
    float const *RESTRICT u_in) {
  auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 3u);
  vec.set(blockIdx, threadIdx);
  u_in += offset;
  if (vec.isValidPosition()) {
    vec.get(0u) = u_in[0u];
    vec.get(1u) = u_in[1u];
    vec.get(2u) = u_in[2u];
  }
}

__global__ void kernel_broadcast(
    gpu::FieldAccessor<float> vec,
    float const *RESTRICT u_in) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    vec.get(0u) = u_in[0u];
    vec.get(1u) = u_in[1u];
    vec.get(2u) = u_in[2u];
  }
}

__global__ void kernel_add(
    gpu::FieldAccessor<float> vec,
    float const *RESTRICT u_in) {
  auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 3u);
  vec.set(blockIdx, threadIdx);
  u_in += offset;
  if (vec.isValidPosition()) {
    vec.get(0u) += u_in[0u];
    vec.get(1u) += u_in[1u];
    vec.get(2u) += u_in[2u];
  }
}

__global__ void kernel_broadcast_add(
    gpu::FieldAccessor<float> vec,
    float const *RESTRICT u_in) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    vec.get(0u) += u_in[0u];
    vec.get(1u) += u_in[1u];
    vec.get(2u) += u_in[2u];
  }
}
// LCOV_EXCL_STOP

Vector3<float> get(
    gpu::GPUField<float> const *vec_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<float> dev_data(3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*vec_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  Vector3<float> vec;
  thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
  return vec;
}

void set(
    gpu::GPUField<float> *vec_field,
    Vector3<float> const &vec,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<float> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const float *>(dev_data_ptr));
  kernel();
}

void add(
    gpu::GPUField<float> *vec_field,
    Vector3<float> const &vec,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<float> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_add);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const float *>(dev_data_ptr));
  kernel();
}

void initialize(
    gpu::GPUField<float> *vec_field,
    Vector3<float> const &vec) {
  CellInterval ci = vec_field->xyzSizeWithGhostLayer();
  thrust::device_vector<float> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_broadcast);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const float *>(dev_data_ptr));
  kernel();
}

void add_to_all(
    gpu::GPUField<float> *vec_field,
    Vector3<float> const &vec) {
  CellInterval ci = vec_field->xyzSizeWithGhostLayer();
  thrust::device_vector<float> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_broadcast_add);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const float *>(dev_data_ptr));
  kernel();
}

std::vector<float> get(
    gpu::GPUField<float> const *vec_field,
    CellInterval const &ci) {
  thrust::device_vector<float> dev_data(ci.numCells() * 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*vec_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<float> out(ci.numCells() * 3u);
  thrust::copy(dev_data.begin(), dev_data.end(), out.data());
  return out;
}

void set(
    gpu::GPUField<float> *vec_field,
    std::vector<float> const &values,
    CellInterval const &ci) {
  thrust::device_vector<float> dev_data(values.begin(), values.end());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const float *>(dev_data_ptr));
  kernel();
}
} // namespace Vector

namespace Flux {
// LCOV_EXCL_START
__global__ void kernel_get(
    gpu::FieldAccessor<float> flux_field,
    float *j_out) {
  auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 13u);
  flux_field.set(blockIdx, threadIdx);
  j_out += offset;
  if (flux_field.isValidPosition()) {
    j_out[0u] = flux_field.get(0u);
    j_out[1u] = flux_field.get(1u);
    j_out[2u] = flux_field.get(2u);
    j_out[3u] = flux_field.get(3u);
    j_out[4u] = flux_field.get(4u);
    j_out[5u] = flux_field.get(5u);
    j_out[6u] = flux_field.get(6u);
    j_out[7u] = flux_field.get(7u);
    j_out[8u] = flux_field.get(8u);
    j_out[9u] = flux_field.get(9u);
    j_out[10u] = flux_field.get(10u);
    j_out[11u] = flux_field.get(11u);
    j_out[12u] = flux_field.get(12u);
  }
}

__global__ void kernel_broadcast(
    gpu::FieldAccessor<float> flux_field,
    float const *RESTRICT j_in) {
  flux_field.set(blockIdx, threadIdx);
  if (flux_field.isValidPosition()) {
    flux_field.get(0u) = j_in[0u];
    flux_field.get(1u) = j_in[1u];
    flux_field.get(2u) = j_in[2u];
    flux_field.get(3u) = j_in[3u];
    flux_field.get(4u) = j_in[4u];
    flux_field.get(5u) = j_in[5u];
    flux_field.get(6u) = j_in[6u];
    flux_field.get(7u) = j_in[7u];
    flux_field.get(8u) = j_in[8u];
    flux_field.get(9u) = j_in[9u];
    flux_field.get(10u) = j_in[10u];
    flux_field.get(11u) = j_in[11u];
    flux_field.get(12u) = j_in[12u];
  }
}

__global__ void kernel_get_vector(
    gpu::FieldAccessor<float> flux_field,
    float *j_out) {
  auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 3u);
  flux_field.set(blockIdx, threadIdx);
  j_out += offset;
  if (flux_field.isValidPosition()) {
    j_out[0u] = float(0.0);
    j_out[1u] = float(0.0);
    j_out[2u] = float(0.0);

    int cx = 0;
    int cy = 0;
    int cz = 0;
    float add_flux;

    cx = 0;
    cy = 1;
    cz = 0;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 1u);
    j_out[1u] += add_flux * 1;
    add_flux = float(0.5) * flux_field.get(1u);
    j_out[1u] += add_flux * -1;
    add_flux = float(0.5) * flux_field.get(0u);
    j_out[0u] += add_flux * -1;
    cx = 1;
    cy = 0;
    cz = 0;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 0u);
    j_out[0u] += add_flux * 1;
    cx = 0;
    cy = 0;
    cz = 1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 2u);
    j_out[2u] += add_flux * 1;
    add_flux = float(0.5) * flux_field.get(2u);
    j_out[2u] += add_flux * -1;
    add_flux = float(0.5) * flux_field.get(4u);
    j_out[0u] += add_flux * -1;
    j_out[1u] += add_flux * 1;
    cx = 1;
    cy = 1;
    cz = 0;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 3u);
    j_out[0u] += add_flux * 1;
    j_out[1u] += add_flux * 1;
    add_flux = float(0.5) * flux_field.get(3u);
    j_out[0u] += add_flux * -1;
    j_out[1u] += add_flux * -1;
    cx = 1;
    cy = -1;
    cz = 0;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 4u);
    j_out[0u] += add_flux * 1;
    j_out[1u] += add_flux * -1;
    cx = 0;
    cy = 1;
    cz = 1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 7u);
    j_out[1u] += add_flux * 1;
    j_out[2u] += add_flux * 1;
    add_flux = float(0.5) * flux_field.get(8u);
    j_out[1u] += add_flux * -1;
    j_out[2u] += add_flux * 1;
    add_flux = float(0.5) * flux_field.get(6u);
    j_out[0u] += add_flux * -1;
    j_out[2u] += add_flux * 1;
    cx = 1;
    cy = 0;
    cz = 1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 5u);
    j_out[0u] += add_flux * 1;
    j_out[2u] += add_flux * 1;
    cx = 0;
    cy = 1;
    cz = -1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 8u);
    j_out[1u] += add_flux * 1;
    j_out[2u] += add_flux * -1;
    add_flux = float(0.5) * flux_field.get(7u);
    j_out[1u] += add_flux * -1;
    j_out[2u] += add_flux * -1;
    add_flux = float(0.5) * flux_field.get(5u);
    j_out[0u] += add_flux * -1;
    j_out[2u] += add_flux * -1;
    cx = 1;
    cy = 0;
    cz = -1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 6u);
    j_out[0u] += add_flux * 1;
    j_out[2u] += add_flux * -1;
    cx = 1;
    cy = 1;
    cz = 1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 9u);
    j_out[0u] += add_flux * 1;
    j_out[1u] += add_flux * 1;
    j_out[2u] += add_flux * 1;
    add_flux = float(0.5) * flux_field.get(12u);
    j_out[0u] += add_flux * -1;
    j_out[1u] += add_flux * 1;
    j_out[2u] += add_flux * 1;
    cx = 1;
    cy = -1;
    cz = 1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 11u);
    j_out[0u] += add_flux * 1;
    j_out[1u] += add_flux * -1;
    j_out[2u] += add_flux * 1;
    add_flux = float(0.5) * flux_field.get(10u);
    j_out[0u] += add_flux * -1;
    j_out[1u] += add_flux * -1;
    j_out[2u] += add_flux * 1;
    cx = 1;
    cy = 1;
    cz = -1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 10u);
    j_out[0u] += add_flux * 1;
    j_out[1u] += add_flux * 1;
    j_out[2u] += add_flux * -1;
    add_flux = float(0.5) * flux_field.get(11u);
    j_out[0u] += add_flux * -1;
    j_out[1u] += add_flux * 1;
    j_out[2u] += add_flux * -1;
    cx = 1;
    cy = -1;
    cz = -1;
    add_flux = float(-0.5) * flux_field.getNeighbor(cx, cy, cz, 12u);
    j_out[0u] += add_flux * 1;
    j_out[1u] += add_flux * -1;
    j_out[2u] += add_flux * -1;
    add_flux = float(0.5) * flux_field.get(9u);
    j_out[0u] += add_flux * -1;
    j_out[1u] += add_flux * -1;
    j_out[2u] += add_flux * -1;
  }
}
// LCOV_EXCL_STOP

std::array<float, 13> get(
    gpu::GPUField<float> const *flux_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<float> dev_data(13u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*flux_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::array<float, 13> vec;
  thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
  return vec;
}

void initialize(
    gpu::GPUField<float> *flux_field,
    std::array<float, 13> const &flux) {
  CellInterval ci = flux_field->xyzSizeWithGhostLayer();
  thrust::device_vector<float> dev_data(flux.data(), flux.data() + 13u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_broadcast);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*flux_field, ci));
  kernel.addParam(const_cast<const float *>(dev_data_ptr));
  kernel();
}

std::vector<float> get(
    gpu::GPUField<float> const *flux_field,
    CellInterval const &ci) {
  thrust::device_vector<float> dev_data(ci.numCells() * 13u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*flux_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<float> out(ci.numCells() * 13u);
  thrust::copy(dev_data.begin(), dev_data.end(), out.data());
  return out;
}

Vector3<float> get_vector(
    gpu::GPUField<float> const *flux_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<float> dev_data(3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get_vector);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*flux_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  Vector3<float> vec;
  thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
  return vec;
}

std::vector<float> get_vector(
    gpu::GPUField<float> const *flux_field,
    CellInterval const &ci) {
  thrust::device_vector<float> dev_data(ci.numCells() * 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get_vector);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<float>::interval(*flux_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<float> out(ci.numCells() * 3u);
  thrust::copy(dev_data.begin(), dev_data.end(), out.data());
  return out;
}
} // namespace Flux

} // namespace accessor
} // namespace ek
} // namespace walberla
