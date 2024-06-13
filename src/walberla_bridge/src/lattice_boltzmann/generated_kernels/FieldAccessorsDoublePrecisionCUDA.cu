/*
 * Copyright (C) 2023-2024 The ESPResSo project
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

// kernel generated with pystencils v1.3.3, lbmpy v1.3.3, lbmpy_walberla/pystencils_walberla from waLBerla commit b0842e1a493ce19ef1bbb8d2cf382fc343970a7f

/**
 * @file
 * Lattice field accessors.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/a16141524c58ab88386e2a0f8fdd7c63c5edd704/python/lbmpy_walberla/templates/LatticeModel.tmpl.h
 */

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Matrix3.h>
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
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // unused variable
#else
#pragma push
#pragma diag_suppress 177 // unused variable
#endif
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#else
// clang compiling CUDA code in host mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#endif
#elif defined(__GNUC__) or defined(__GNUG__)
#define RESTRICT __restrict__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
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
namespace lbm {
namespace accessor {

namespace Population {
__global__ void kernel_get_interval(
    gpu::FieldAccessor<double> pdf,
    double *RESTRICT const pop) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 19u);
    pop[offset + 0u] = pdf.get(0u);
    pop[offset + 1u] = pdf.get(1u);
    pop[offset + 2u] = pdf.get(2u);
    pop[offset + 3u] = pdf.get(3u);
    pop[offset + 4u] = pdf.get(4u);
    pop[offset + 5u] = pdf.get(5u);
    pop[offset + 6u] = pdf.get(6u);
    pop[offset + 7u] = pdf.get(7u);
    pop[offset + 8u] = pdf.get(8u);
    pop[offset + 9u] = pdf.get(9u);
    pop[offset + 10u] = pdf.get(10u);
    pop[offset + 11u] = pdf.get(11u);
    pop[offset + 12u] = pdf.get(12u);
    pop[offset + 13u] = pdf.get(13u);
    pop[offset + 14u] = pdf.get(14u);
    pop[offset + 15u] = pdf.get(15u);
    pop[offset + 16u] = pdf.get(16u);
    pop[offset + 17u] = pdf.get(17u);
    pop[offset + 18u] = pdf.get(18u);
  }
}

__global__ void kernel_get(
    gpu::FieldAccessor<double> pdf,
    double *RESTRICT const pop) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    pop[0u] = pdf.get(0u);
    pop[1u] = pdf.get(1u);
    pop[2u] = pdf.get(2u);
    pop[3u] = pdf.get(3u);
    pop[4u] = pdf.get(4u);
    pop[5u] = pdf.get(5u);
    pop[6u] = pdf.get(6u);
    pop[7u] = pdf.get(7u);
    pop[8u] = pdf.get(8u);
    pop[9u] = pdf.get(9u);
    pop[10u] = pdf.get(10u);
    pop[11u] = pdf.get(11u);
    pop[12u] = pdf.get(12u);
    pop[13u] = pdf.get(13u);
    pop[14u] = pdf.get(14u);
    pop[15u] = pdf.get(15u);
    pop[16u] = pdf.get(16u);
    pop[17u] = pdf.get(17u);
    pop[18u] = pdf.get(18u);
  }
}

__global__ void kernel_set_interval(
    gpu::FieldAccessor<double> pdf,
    double const *RESTRICT const pop) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 19u);
    pdf.get(0u) = pop[offset + 0u];
    pdf.get(1u) = pop[offset + 1u];
    pdf.get(2u) = pop[offset + 2u];
    pdf.get(3u) = pop[offset + 3u];
    pdf.get(4u) = pop[offset + 4u];
    pdf.get(5u) = pop[offset + 5u];
    pdf.get(6u) = pop[offset + 6u];
    pdf.get(7u) = pop[offset + 7u];
    pdf.get(8u) = pop[offset + 8u];
    pdf.get(9u) = pop[offset + 9u];
    pdf.get(10u) = pop[offset + 10u];
    pdf.get(11u) = pop[offset + 11u];
    pdf.get(12u) = pop[offset + 12u];
    pdf.get(13u) = pop[offset + 13u];
    pdf.get(14u) = pop[offset + 14u];
    pdf.get(15u) = pop[offset + 15u];
    pdf.get(16u) = pop[offset + 16u];
    pdf.get(17u) = pop[offset + 17u];
    pdf.get(18u) = pop[offset + 18u];
  }
}

__global__ void kernel_set(
    gpu::FieldAccessor<double> pdf,
    double const *RESTRICT const pop) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    pdf.get(0u) = pop[0u];
    pdf.get(1u) = pop[1u];
    pdf.get(2u) = pop[2u];
    pdf.get(3u) = pop[3u];
    pdf.get(4u) = pop[4u];
    pdf.get(5u) = pop[5u];
    pdf.get(6u) = pop[6u];
    pdf.get(7u) = pop[7u];
    pdf.get(8u) = pop[8u];
    pdf.get(9u) = pop[9u];
    pdf.get(10u) = pop[10u];
    pdf.get(11u) = pop[11u];
    pdf.get(12u) = pop[12u];
    pdf.get(13u) = pop[13u];
    pdf.get(14u) = pop[14u];
    pdf.get(15u) = pop[15u];
    pdf.get(16u) = pop[16u];
    pdf.get(17u) = pop[17u];
    pdf.get(18u) = pop[18u];
  }
}

std::array<double, 19u> get(
    gpu::GPUField<double> const *pdf_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(19u, double{0});
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::array<double, 19u> pop;
  thrust::copy(dev_data.begin(), dev_data.end(), pop.data());
  return pop;
}

void set(
    gpu::GPUField<double> *pdf_field,
    std::array<double, 19u> const &pop,
    Cell const &cell) {
  thrust::device_vector<double> dev_data(pop.data(), pop.data() + 19u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  CellInterval ci(cell, cell);
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

void initialize(
    gpu::GPUField<double> *pdf_field,
    std::array<double, 19u> const &pop) {
  CellInterval ci = pdf_field->xyzSizeWithGhostLayer();
  thrust::device_vector<double> dev_data(pop.data(), pop.data() + 19u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

std::vector<double> get(
    gpu::GPUField<double> const *pdf_field,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(ci.numCells() * 19u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get_interval);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<double> out(ci.numCells() * 19u);
  thrust::copy(dev_data.begin(), dev_data.end(), out.data());
  return out;
}

void set(
    gpu::GPUField<double> *pdf_field,
    std::vector<double> const &values,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(values.begin(), values.end());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set_interval);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}
} // namespace Population

namespace Vector {
__global__ void kernel_get_interval(
    gpu::FieldAccessor<double> vec,
    double *const out) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 3u);
    out[offset + 0u] = vec.get(0u);
    out[offset + 1u] = vec.get(1u);
    out[offset + 2u] = vec.get(2u);
  }
}

__global__ void kernel_get(
    gpu::FieldAccessor<double> vec,
    double *const out) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    out[0u] = vec.get(0u);
    out[1u] = vec.get(1u);
    out[2u] = vec.get(2u);
  }
}

__global__ void kernel_set_interval(
    gpu::FieldAccessor<double> vec,
    double const *RESTRICT const u) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 3u);
    vec.get(0u) = u[offset + 0u];
    vec.get(1u) = u[offset + 1u];
    vec.get(2u) = u[offset + 2u];
  }
}

__global__ void kernel_set(
    gpu::FieldAccessor<double> vec,
    const double *RESTRICT const u) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    vec.get(0u) = u[0u];
    vec.get(1u) = u[1u];
    vec.get(2u) = u[2u];
  }
}

__global__ void kernel_add_interval(
    gpu::FieldAccessor<double> vec,
    double const *RESTRICT const u) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 3u);
    vec.get(0u) += u[offset + 0u];
    vec.get(1u) += u[offset + 1u];
    vec.get(2u) += u[offset + 2u];
  }
}

__global__ void kernel_add(
    gpu::FieldAccessor<double> vec,
    double const *RESTRICT const u) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    vec.get(0u) += u[0u];
    vec.get(1u) += u[1u];
    vec.get(2u) += u[2u];
  }
}

Vector3<double> get(
    gpu::GPUField<double> const *vec_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  Vector3<double> vec;
  thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
  return vec;
}

void set(
    gpu::GPUField<double> *vec_field,
    Vector3<double> const &vec,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

void add(
    gpu::GPUField<double> *vec_field,
    Vector3<double> const &vec,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_add);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

void initialize(
    gpu::GPUField<double> *vec_field,
    Vector3<double> const &vec) {
  CellInterval ci = vec_field->xyzSizeWithGhostLayer();
  thrust::device_vector<double> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

void add_to_all(
    gpu::GPUField<double> *vec_field,
    Vector3<double> const &vec) {
  CellInterval ci = vec_field->xyzSizeWithGhostLayer();
  thrust::device_vector<double> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_add);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

std::vector<double> get(
    gpu::GPUField<double> const *vec_field,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(ci.numCells() * 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get_interval);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<double> out(ci.numCells() * 3u);
  thrust::copy(dev_data.begin(), dev_data.end(), out.data());
  return out;
}

void set(
    gpu::GPUField<double> *vec_field,
    std::vector<double> const &values,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(values.begin(), values.end());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set_interval);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}
} // namespace Vector

namespace Interpolation {
/** @brief Calculate interpolation weights. */
static __forceinline__ __device__ void calculate_weights(
    double const *RESTRICT const pos,
    int *RESTRICT const corner,
    double *RESTRICT const weights,
    uint gl) {
#pragma unroll
  for (int dim = 0; dim < 3; ++dim) {
    auto const fractional_index = pos[dim] - double{0.5};
    auto const nmp = floorf(fractional_index);
    auto const distance = fractional_index - nmp - double{0.5};
    corner[dim] = __double2int_rn(nmp) + static_cast<int>(gl);
    weights[dim * 2 + 0] = double{0.5} - distance;
    weights[dim * 2 + 1] = double{0.5} + distance;
  }
}

__global__ void kernel_get(
    gpu::FieldAccessor<double> vec,
    double const *RESTRICT const pos,
    double *RESTRICT const vel,
    uint n_pos,
    uint gl) {

  uint pos_index = blockIdx.y * gridDim.x * blockDim.x +
                   blockDim.x * blockIdx.x + threadIdx.x;

  vec.set({0u, 0u, 0u}, {0u, 0u, 0u});
  if (vec.isValidPosition() and pos_index < n_pos) {
    auto const array_offset = pos_index * uint(3u);
    int corner[3];
    double weights[3][2];
    calculate_weights(pos + array_offset, corner, &weights[0][0], gl);
#pragma unroll
    for (int i = 0; i < 2; i++) {
      auto const cx = corner[0] + i;
      auto const wx = weights[0][i];
#pragma unroll
      for (int j = 0; j < 2; j++) {
        auto const cy = corner[1] + j;
        auto const wxy = wx * weights[1][j];
#pragma unroll
        for (int k = 0; k < 2; k++) {
          auto const cz = corner[2] + k;
          auto const weight = wxy * weights[2][k];
          vel[array_offset + 0u] += weight * vec.getNeighbor(cx, cy, cz, 0u);
          vel[array_offset + 1u] += weight * vec.getNeighbor(cx, cy, cz, 1u);
          vel[array_offset + 2u] += weight * vec.getNeighbor(cx, cy, cz, 2u);
        }
      }
    }
  }
}

__global__ void kernel_set(
    gpu::FieldAccessor<double> vec,
    double const *RESTRICT const pos,
    double const *RESTRICT const forces,
    uint n_pos,
    uint gl) {

  uint pos_index = blockIdx.y * gridDim.x * blockDim.x +
                   blockDim.x * blockIdx.x + threadIdx.x;

  vec.set({0u, 0u, 0u}, {0u, 0u, 0u});
  if (vec.isValidPosition() and pos_index < n_pos) {
    auto const array_offset = pos_index * uint(3u);
    int corner[3];
    double weights[3][2];
    calculate_weights(pos + array_offset, corner, &weights[0][0], gl);
#pragma unroll
    for (int i = 0; i < 2; i++) {
      auto const cx = corner[0] + i;
      auto const wx = weights[0][i];
#pragma unroll
      for (int j = 0; j < 2; j++) {
        auto const cy = corner[1] + j;
        auto const wxy = wx * weights[1][j];
#pragma unroll
        for (int k = 0; k < 2; k++) {
          auto const cz = corner[2] + k;
          auto const weight = wxy * weights[2][k];
          atomicAdd(&vec.getNeighbor(cx, cy, cz, 0u),
                    weight * forces[array_offset + 0u]);
          atomicAdd(&vec.getNeighbor(cx, cy, cz, 1u),
                    weight * forces[array_offset + 1u]);
          atomicAdd(&vec.getNeighbor(cx, cy, cz, 2u),
                    weight * forces[array_offset + 2u]);
        }
      }
    }
  }
}

static dim3 calculate_dim_grid(uint const threads_x,
                               uint const blocks_per_grid_y,
                               uint const threads_per_block) {
  assert(threads_x >= 1u);
  assert(blocks_per_grid_y >= 1u);
  assert(threads_per_block >= 1u);
  auto const threads_y = threads_per_block * blocks_per_grid_y;
  auto const blocks_per_grid_x = (threads_x + threads_y - 1) / threads_y;
  return make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
}

std::vector<double>
get(
    gpu::GPUField<double> const *vec_field,
    std::vector<double> const &pos,
    uint gl) {
  thrust::device_vector<double> dev_pos(pos.begin(), pos.end());
  thrust::device_vector<double> dev_vel(pos.size());
  auto const dev_pos_ptr = thrust::raw_pointer_cast(dev_pos.data());
  auto const dev_vel_ptr = thrust::raw_pointer_cast(dev_vel.data());

  auto const threads_per_block = uint(64u);
  auto const n_pos = static_cast<uint>(pos.size() / 3ul);
  auto const dim_grid = calculate_dim_grid(n_pos, 4u, threads_per_block);
  kernel_get<<<dim_grid, threads_per_block, 0u, nullptr>>>(
      gpu::FieldIndexing<double>::withGhostLayerXYZ(*vec_field, gl).gpuAccess(),
      dev_pos_ptr, dev_vel_ptr, n_pos, gl);

  std::vector<double> out(pos.size());
  thrust::copy(dev_vel.begin(), dev_vel.end(), out.data());
  return out;
}

void set(
    gpu::GPUField<double> const *vec_field,
    std::vector<double> const &pos,
    std::vector<double> const &forces,
    uint gl) {
  thrust::device_vector<double> dev_pos(pos.begin(), pos.end());
  thrust::device_vector<double> dev_for(forces.begin(), forces.end());
  auto const dev_pos_ptr = thrust::raw_pointer_cast(dev_pos.data());
  auto const dev_for_ptr = thrust::raw_pointer_cast(dev_for.data());

  auto const threads_per_block = uint(64u);
  auto const n_pos = static_cast<uint>(pos.size() / 3ul);
  auto const dim_grid = calculate_dim_grid(n_pos, 4u, threads_per_block);
  kernel_set<<<dim_grid, threads_per_block, 0u, nullptr>>>(
      gpu::FieldIndexing<double>::withGhostLayerXYZ(*vec_field, gl).gpuAccess(),
      dev_pos_ptr, dev_for_ptr, n_pos, gl);
}
} // namespace Interpolation

namespace Equilibrium {
__device__ void kernel_set_device(
    gpu::FieldAccessor<double> pdf,
    double const *RESTRICT const u,
    double rho) {

  pdf.get(0u) = rho * -0.33333333333333331 * (u[0] * u[0]) + rho * -0.33333333333333331 * (u[1] * u[1]) + rho * -0.33333333333333331 * (u[2] * u[2]) + rho * 0.33333333333333331;
  pdf.get(1u) = rho * -0.16666666666666666 * (u[0] * u[0]) + rho * -0.16666666666666666 * (u[2] * u[2]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[1] + rho * 0.16666666666666666 * (u[1] * u[1]);
  pdf.get(2u) = rho * -0.16666666666666666 * u[1] + rho * -0.16666666666666666 * (u[0] * u[0]) + rho * -0.16666666666666666 * (u[2] * u[2]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u[1] * u[1]);
  pdf.get(3u) = rho * -0.16666666666666666 * u[0] + rho * -0.16666666666666666 * (u[1] * u[1]) + rho * -0.16666666666666666 * (u[2] * u[2]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u[0] * u[0]);
  pdf.get(4u) = rho * -0.16666666666666666 * (u[1] * u[1]) + rho * -0.16666666666666666 * (u[2] * u[2]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[0] + rho * 0.16666666666666666 * (u[0] * u[0]);
  pdf.get(5u) = rho * -0.16666666666666666 * (u[0] * u[0]) + rho * -0.16666666666666666 * (u[1] * u[1]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[2] + rho * 0.16666666666666666 * (u[2] * u[2]);
  pdf.get(6u) = rho * -0.16666666666666666 * u[2] + rho * -0.16666666666666666 * (u[0] * u[0]) + rho * -0.16666666666666666 * (u[1] * u[1]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u[2] * u[2]);
  pdf.get(7u) = rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[1] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[1] * u[1]);
  pdf.get(8u) = rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] + rho * 0.083333333333333329 * u[1] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.25 * u[0] * u[1];
  pdf.get(9u) = rho * -0.083333333333333329 * u[0] + rho * -0.083333333333333329 * u[1] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.25 * u[0] * u[1];
  pdf.get(10u) = rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[0] * u[1] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[1] * u[1]);
  pdf.get(11u) = rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] + rho * 0.083333333333333329 * u[2] + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[1] * u[2];
  pdf.get(12u) = rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[1] * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[2] + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf.get(13u) = rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[2] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf.get(14u) = rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] + rho * 0.083333333333333329 * u[2] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[0] * u[2];
  pdf.get(15u) = rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[1] * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf.get(16u) = rho * -0.083333333333333329 * u[1] + rho * -0.083333333333333329 * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[1] * u[2];
  pdf.get(17u) = rho * -0.083333333333333329 * u[0] + rho * -0.083333333333333329 * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[0] * u[2];
  pdf.get(18u) = rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[0] * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[2] * u[2]);
}
} // namespace Equilibrium

namespace Density {
__global__ void kernel_get(
    gpu::FieldAccessor<double> pdf,
    double *RESTRICT const out) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, uint(1u));
    double const f_0 = pdf.get(0u);
    double const f_1 = pdf.get(1u);
    double const f_2 = pdf.get(2u);
    double const f_3 = pdf.get(3u);
    double const f_4 = pdf.get(4u);
    double const f_5 = pdf.get(5u);
    double const f_6 = pdf.get(6u);
    double const f_7 = pdf.get(7u);
    double const f_8 = pdf.get(8u);
    double const f_9 = pdf.get(9u);
    double const f_10 = pdf.get(10u);
    double const f_11 = pdf.get(11u);
    double const f_12 = pdf.get(12u);
    double const f_13 = pdf.get(13u);
    double const f_14 = pdf.get(14u);
    double const f_15 = pdf.get(15u);
    double const f_16 = pdf.get(16u);
    double const f_17 = pdf.get(17u);
    double const f_18 = pdf.get(18u);
    const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const double vel1Term = f_1 + f_11 + f_15 + f_7;
    const double vel2Term = f_12 + f_13 + f_5;
    const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
    out[offset] = rho;
  }
}

__global__ void kernel_set(
    gpu::FieldAccessor<double> pdf,
    double const *RESTRICT const rho_in) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, uint(1u));
    double const f_0 = pdf.get(0u);
    double const f_1 = pdf.get(1u);
    double const f_2 = pdf.get(2u);
    double const f_3 = pdf.get(3u);
    double const f_4 = pdf.get(4u);
    double const f_5 = pdf.get(5u);
    double const f_6 = pdf.get(6u);
    double const f_7 = pdf.get(7u);
    double const f_8 = pdf.get(8u);
    double const f_9 = pdf.get(9u);
    double const f_10 = pdf.get(10u);
    double const f_11 = pdf.get(11u);
    double const f_12 = pdf.get(12u);
    double const f_13 = pdf.get(13u);
    double const f_14 = pdf.get(14u);
    double const f_15 = pdf.get(15u);
    double const f_16 = pdf.get(16u);
    double const f_17 = pdf.get(17u);
    double const f_18 = pdf.get(18u);
    const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
    const double vel1Term = f_1 + f_11 + f_15 + f_7;
    const double momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
    const double vel2Term = f_12 + f_13 + f_5;
    const double momdensity_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
    const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;

    // calculate current velocity (before density change)
    double const conversion = double(1) / rho;
    double const u_old[3] = {momdensity_0 * conversion, momdensity_1 * conversion, momdensity_2 * conversion};

    Equilibrium::kernel_set_device(pdf, u_old, rho_in[offset]);
  }
}

double get(
    gpu::GPUField<double> const *pdf_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(1u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  double rho = dev_data[0u];
  return rho;
}

void set(
    gpu::GPUField<double> *pdf_field,
    const double rho,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(1u, rho);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

std::vector<double> get(
    gpu::GPUField<double> const *pdf_field,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(ci.numCells());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<double> out(ci.numCells());
  thrust::copy(dev_data.begin(), dev_data.end(), out.begin());
  return out;
}

void set(
    gpu::GPUField<double> *pdf_field,
    std::vector<double> const &values,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(values.begin(), values.end());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}
} // namespace Density

namespace Velocity {
__global__ void kernel_set(
    gpu::FieldAccessor<double> pdf,
    gpu::FieldAccessor<double> force,
    double const *RESTRICT const u_in) {
  pdf.set(blockIdx, threadIdx);
  force.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, uint(3u));
    uint const bufsize = 3u;
    double const *RESTRICT const u = u_in + bufsize * offset;
    double const f_0 = pdf.get(0u);
    double const f_1 = pdf.get(1u);
    double const f_2 = pdf.get(2u);
    double const f_3 = pdf.get(3u);
    double const f_4 = pdf.get(4u);
    double const f_5 = pdf.get(5u);
    double const f_6 = pdf.get(6u);
    double const f_7 = pdf.get(7u);
    double const f_8 = pdf.get(8u);
    double const f_9 = pdf.get(9u);
    double const f_10 = pdf.get(10u);
    double const f_11 = pdf.get(11u);
    double const f_12 = pdf.get(12u);
    double const f_13 = pdf.get(13u);
    double const f_14 = pdf.get(14u);
    double const f_15 = pdf.get(15u);
    double const f_16 = pdf.get(16u);
    double const f_17 = pdf.get(17u);
    double const f_18 = pdf.get(18u);
    const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const double vel1Term = f_1 + f_11 + f_15 + f_7;
    const double vel2Term = f_12 + f_13 + f_5;
    const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
    const double u_0 = -force.get(0) * 0.50000000000000000 / rho + u[0];
    const double u_1 = -force.get(1) * 0.50000000000000000 / rho + u[1];
    const double u_2 = -force.get(2) * 0.50000000000000000 / rho + u[2];
    double u_new[3] = {u_0, u_1, u_2};

    Equilibrium::kernel_set_device(pdf, u_new, rho);
  }
}

void set(
    gpu::GPUField<double> *pdf_field,
    gpu::GPUField<double> *force_field,
    Vector3<double> const &u,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(u.data(), u.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*force_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}
} // namespace Velocity

namespace MomentumDensity {
__global__ void kernel_sum(
    gpu::FieldAccessor<double> pdf,
    gpu::FieldAccessor<double> force,
    double *RESTRICT const out) {
  pdf.set(blockIdx, threadIdx);
  force.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    uint const bufsize = 3u;
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, bufsize);
    double const f_0 = pdf.get(0u);
    double const f_1 = pdf.get(1u);
    double const f_2 = pdf.get(2u);
    double const f_3 = pdf.get(3u);
    double const f_4 = pdf.get(4u);
    double const f_5 = pdf.get(5u);
    double const f_6 = pdf.get(6u);
    double const f_7 = pdf.get(7u);
    double const f_8 = pdf.get(8u);
    double const f_9 = pdf.get(9u);
    double const f_10 = pdf.get(10u);
    double const f_11 = pdf.get(11u);
    double const f_12 = pdf.get(12u);
    double const f_13 = pdf.get(13u);
    double const f_14 = pdf.get(14u);
    double const f_15 = pdf.get(15u);
    double const f_16 = pdf.get(16u);
    double const f_17 = pdf.get(17u);
    double const f_18 = pdf.get(18u);
    const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
    const double vel1Term = f_1 + f_11 + f_15 + f_7;
    const double momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
    const double vel2Term = f_12 + f_13 + f_5;
    const double momdensity_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
    const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
    const double md_0 = force.get(0) * 0.50000000000000000 + momdensity_0;
    const double md_1 = force.get(1) * 0.50000000000000000 + momdensity_1;
    const double md_2 = force.get(2) * 0.50000000000000000 + momdensity_2;
    out[bufsize * offset + 0u] += md_0;
    out[bufsize * offset + 1u] += md_1;
    out[bufsize * offset + 2u] += md_2;
  }
}

Vector3<double> reduce(
    gpu::GPUField<double> const *pdf_field,
    gpu::GPUField<double> const *force_field) {
  thrust::device_vector<double> dev_data(3u, double{0});
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
    Cell cell(x, y, z);
    CellInterval ci(cell, cell);
    auto kernel = gpu::make_kernel(kernel_sum);
    kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
    kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*force_field, ci));
    kernel.addParam(dev_data_ptr);
    kernel();
  });
  Vector3<double> mom(double{0});
  thrust::copy(dev_data.begin(), dev_data.begin() + 3u, mom.data());
  return mom;
}
} // namespace MomentumDensity

namespace PressureTensor {
__global__ void kernel_get(
    gpu::FieldAccessor<double> pdf,
    double *RESTRICT const out) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    uint const bufsize = 9u;
    uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, bufsize);
    double const f_0 = pdf.get(0u);
    double const f_1 = pdf.get(1u);
    double const f_2 = pdf.get(2u);
    double const f_3 = pdf.get(3u);
    double const f_4 = pdf.get(4u);
    double const f_5 = pdf.get(5u);
    double const f_6 = pdf.get(6u);
    double const f_7 = pdf.get(7u);
    double const f_8 = pdf.get(8u);
    double const f_9 = pdf.get(9u);
    double const f_10 = pdf.get(10u);
    double const f_11 = pdf.get(11u);
    double const f_12 = pdf.get(12u);
    double const f_13 = pdf.get(13u);
    double const f_14 = pdf.get(14u);
    double const f_15 = pdf.get(15u);
    double const f_16 = pdf.get(16u);
    double const f_17 = pdf.get(17u);
    double const f_18 = pdf.get(18u);
    const double p_0 = f_10 + f_13 + f_14 + f_17 + f_18 + f_3 + f_4 + f_7 + f_8 + f_9;
    const double p_1 = -f_10 - f_7 + f_8 + f_9;
    const double p_2 = -f_13 + f_14 + f_17 - f_18;
    const double p_3 = -f_10 - f_7 + f_8 + f_9;
    const double p_4 = f_1 + f_10 + f_11 + f_12 + f_15 + f_16 + f_2 + f_7 + f_8 + f_9;
    const double p_5 = f_11 - f_12 - f_15 + f_16;
    const double p_6 = -f_13 + f_14 + f_17 - f_18;
    const double p_7 = f_11 - f_12 - f_15 + f_16;
    const double p_8 = f_11 + f_12 + f_13 + f_14 + f_15 + f_16 + f_17 + f_18 + f_5 + f_6;
    out[bufsize * offset + 0u] = p_0;
    out[bufsize * offset + 1u] = p_1;
    out[bufsize * offset + 2u] = p_2;

    out[bufsize * offset + 3u] = p_3;
    out[bufsize * offset + 4u] = p_4;
    out[bufsize * offset + 5u] = p_5;

    out[bufsize * offset + 6u] = p_6;
    out[bufsize * offset + 7u] = p_7;
    out[bufsize * offset + 8u] = p_8;
  }
}

Matrix3<double> get(
    gpu::GPUField<double> const *pdf_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(9u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = gpu::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  Matrix3<double> out;
  thrust::copy(dev_data.begin(), dev_data.begin() + 9u, out.data());
  return out;
}
} // namespace PressureTensor

} // namespace accessor
} // namespace lbm
} // namespace walberla
