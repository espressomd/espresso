/*
 * Copyright (C) 2023 The ESPResSo project
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

// kernel generated with pystencils v1.2, lbmpy v1.2, lbmpy_walberla/pystencils_walberla from waLBerla commit 065ce5f311850371a97ac4766f47dbb5ca8424ba

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

#include <cuda/FieldAccessor.h>
#include <cuda/FieldIndexing.h>
#include <cuda/GPUField.h>
#include <cuda/Kernel.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

#include <array>
#include <tuple>
#include <vector>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

__device__ inline uint get_num_threads(uint3 gridDim, uint3 blockDim) {
  return gridDim.x * gridDim.y * gridDim.z * blockDim.x * blockDim.y * blockDim.z;
}

__device__ inline uint getLinearIndexXYZF(uint3 blockIdx, uint3 threadIdx, uint3 gridDim, uint3 blockDim) {
  auto const x = threadIdx.x;
  auto const y = blockIdx.x;
  auto const z = blockIdx.y;
  auto const f = blockIdx.z;
  auto const xSize = blockDim.x;
  auto const ySize = gridDim.x;
  auto const zSize = gridDim.y;
  return x +
         y * xSize +
         z * xSize * ySize +
         f * xSize * ySize * zSize;
}

__device__ inline uint getLinearIndexFZYX(uint3 blockIdx, uint3 threadIdx, uint3 gridDim, uint3 blockDim, uint fOffset) {
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
__global__ void kernel_get(
    cuda::FieldAccessor<double> pdf,
    double *RESTRICT const pop) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, 19u);
    pop[offset + 0u] = pdf.get(0);
    pop[offset + 1u] = pdf.get(1);
    pop[offset + 2u] = pdf.get(2);
    pop[offset + 3u] = pdf.get(3);
    pop[offset + 4u] = pdf.get(4);
    pop[offset + 5u] = pdf.get(5);
    pop[offset + 6u] = pdf.get(6);
    pop[offset + 7u] = pdf.get(7);
    pop[offset + 8u] = pdf.get(8);
    pop[offset + 9u] = pdf.get(9);
    pop[offset + 10u] = pdf.get(10);
    pop[offset + 11u] = pdf.get(11);
    pop[offset + 12u] = pdf.get(12);
    pop[offset + 13u] = pdf.get(13);
    pop[offset + 14u] = pdf.get(14);
    pop[offset + 15u] = pdf.get(15);
    pop[offset + 16u] = pdf.get(16);
    pop[offset + 17u] = pdf.get(17);
    pop[offset + 18u] = pdf.get(18);
  }
}

__global__ void kernel_set(
    cuda::FieldAccessor<double> pdf,
    const double *RESTRICT const pop) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, 19u);
    pdf.get(0) = pop[offset + 0u];
    pdf.get(1) = pop[offset + 1u];
    pdf.get(2) = pop[offset + 2u];
    pdf.get(3) = pop[offset + 3u];
    pdf.get(4) = pop[offset + 4u];
    pdf.get(5) = pop[offset + 5u];
    pdf.get(6) = pop[offset + 6u];
    pdf.get(7) = pop[offset + 7u];
    pdf.get(8) = pop[offset + 8u];
    pdf.get(9) = pop[offset + 9u];
    pdf.get(10) = pop[offset + 10u];
    pdf.get(11) = pop[offset + 11u];
    pdf.get(12) = pop[offset + 12u];
    pdf.get(13) = pop[offset + 13u];
    pdf.get(14) = pop[offset + 14u];
    pdf.get(15) = pop[offset + 15u];
    pdf.get(16) = pop[offset + 16u];
    pdf.get(17) = pop[offset + 17u];
    pdf.get(18) = pop[offset + 18u];
  }
}

std::array<double, 19u> get(
    cuda::GPUField<double> const *pdf_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(19u, double{0});
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::array<double, 19u> pop;
  thrust::copy(dev_data.begin(), dev_data.end(), pop.data());
  return pop;
}

void set(
    cuda::GPUField<double> *pdf_field,
    std::array<double, 19u> const &pop,
    Cell const &cell) {
  thrust::device_vector<double> dev_data(pop.data(), pop.data() + 19u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  CellInterval ci(cell, cell);
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

void broadcast(
    cuda::GPUField<double> *pdf_field,
    std::array<double, 19u> const &pop) {
  CellInterval ci = pdf_field->xyzSizeWithGhostLayer();
  thrust::device_vector<double> dev_data(pop.data(), pop.data() + 19u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

std::vector<double> get(
    cuda::GPUField<double> const *pdf_field,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(ci.numCells() * 19u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<double> out(ci.numCells() * 19u);
  thrust::copy(dev_data.begin(), dev_data.end(), out.data());
  return out;
}

void set(
    cuda::GPUField<double> *pdf_field,
    std::vector<double> const &values,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(values.begin(), values.end());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}
} // namespace Population

namespace Vector {
__global__ void kernel_get(
    cuda::FieldAccessor<double> vec,
    double *const out) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, 3u);
    out[offset + 0u] = vec.get(0);
    out[offset + 1u] = vec.get(1);
    out[offset + 2u] = vec.get(2);
  }
}

__global__ void kernel_set(
    cuda::FieldAccessor<double> vec,
    const double *RESTRICT const u) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, 3u);
    vec.get(0) = u[offset + 0u];
    vec.get(1) = u[offset + 1u];
    vec.get(2) = u[offset + 2u];
  }
}

__global__ void kernel_add(
    cuda::FieldAccessor<double> vec,
    const double *RESTRICT const u) {
  vec.set(blockIdx, threadIdx);
  if (vec.isValidPosition()) {
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, 3u);
    vec.get(0) += u[offset + 0u];
    vec.get(1) += u[offset + 1u];
    vec.get(2) += u[offset + 2u];
  }
}

Vector3<double> get(
    cuda::GPUField<double> const *vec_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  Vector3<double> vec;
  thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
  return vec;
}

void set(
    cuda::GPUField<double> *vec_field,
    Vector3<double> const &vec,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

void add(
    cuda::GPUField<double> *vec_field,
    Vector3<double> const &vec,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_add);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

void broadcast(
    cuda::GPUField<double> *vec_field,
    Vector3<double> const &vec) {
  CellInterval ci = vec_field->xyzSizeWithGhostLayer();
  thrust::device_vector<double> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

void add_to_all(
    cuda::GPUField<double> *vec_field,
    Vector3<double> const &vec) {
  CellInterval ci = vec_field->xyzSizeWithGhostLayer();
  thrust::device_vector<double> dev_data(vec.data(), vec.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_add);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

std::vector<double> get(
    cuda::GPUField<double> const *vec_field,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(ci.numCells() * 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<double> out(ci.numCells() * 3u);
  thrust::copy(dev_data.begin(), dev_data.end(), out.data());
  return out;
}

void set(
    cuda::GPUField<double> *vec_field,
    std::vector<double> const &values,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(values.begin(), values.end());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*vec_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}
} // namespace Vector

namespace Equilibrium {
__device__ void kernel_set_device(
    cuda::FieldAccessor<double> pdf,
    const double *RESTRICT const u,
    double rho) {

  pdf.get(0) = rho * -0.33333333333333331 * (u[0] * u[0]) + rho * -0.33333333333333331 * (u[1] * u[1]) + rho * -0.33333333333333331 * (u[2] * u[2]) + rho * 0.33333333333333331;
  pdf.get(1) = rho * -0.16666666666666666 * (u[0] * u[0]) + rho * -0.16666666666666666 * (u[2] * u[2]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[1] + rho * 0.16666666666666666 * (u[1] * u[1]);
  pdf.get(2) = rho * -0.16666666666666666 * u[1] + rho * -0.16666666666666666 * (u[0] * u[0]) + rho * -0.16666666666666666 * (u[2] * u[2]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u[1] * u[1]);
  pdf.get(3) = rho * -0.16666666666666666 * u[0] + rho * -0.16666666666666666 * (u[1] * u[1]) + rho * -0.16666666666666666 * (u[2] * u[2]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u[0] * u[0]);
  pdf.get(4) = rho * -0.16666666666666666 * (u[1] * u[1]) + rho * -0.16666666666666666 * (u[2] * u[2]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[0] + rho * 0.16666666666666666 * (u[0] * u[0]);
  pdf.get(5) = rho * -0.16666666666666666 * (u[0] * u[0]) + rho * -0.16666666666666666 * (u[1] * u[1]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * u[2] + rho * 0.16666666666666666 * (u[2] * u[2]);
  pdf.get(6) = rho * -0.16666666666666666 * u[2] + rho * -0.16666666666666666 * (u[0] * u[0]) + rho * -0.16666666666666666 * (u[1] * u[1]) + rho * 0.055555555555555552 + rho * 0.16666666666666666 * (u[2] * u[2]);
  pdf.get(7) = rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[1] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[1] * u[1]);
  pdf.get(8) = rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] + rho * 0.083333333333333329 * u[1] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.25 * u[0] * u[1];
  pdf.get(9) = rho * -0.083333333333333329 * u[0] + rho * -0.083333333333333329 * u[1] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.25 * u[0] * u[1];
  pdf.get(10) = rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[0] * u[1] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[1] * u[1]);
  pdf.get(11) = rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] + rho * 0.083333333333333329 * u[2] + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[1] * u[2];
  pdf.get(12) = rho * -0.083333333333333329 * u[1] + rho * -0.25 * u[1] * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[2] + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf.get(13) = rho * -0.083333333333333329 * u[0] + rho * -0.25 * u[0] * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[2] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf.get(14) = rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] + rho * 0.083333333333333329 * u[2] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[0] * u[2];
  pdf.get(15) = rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[1] * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[1] + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.083333333333333329 * (u[2] * u[2]);
  pdf.get(16) = rho * -0.083333333333333329 * u[1] + rho * -0.083333333333333329 * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[1] * u[1]) + rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[1] * u[2];
  pdf.get(17) = rho * -0.083333333333333329 * u[0] + rho * -0.083333333333333329 * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[2] * u[2]) + rho * 0.25 * u[0] * u[2];
  pdf.get(18) = rho * -0.083333333333333329 * u[2] + rho * -0.25 * u[0] * u[2] + rho * 0.027777777777777776 + rho * 0.083333333333333329 * u[0] + rho * 0.083333333333333329 * (u[0] * u[0]) + rho * 0.083333333333333329 * (u[2] * u[2]);
}
} // namespace Equilibrium

namespace Density {
__global__ void kernel_get(
    cuda::FieldAccessor<double> pdf,
    double *RESTRICT const out) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, uint(1u));
    const double f_0 = pdf.get(0);
    const double f_1 = pdf.get(1);
    const double f_2 = pdf.get(2);
    const double f_3 = pdf.get(3);
    const double f_4 = pdf.get(4);
    const double f_5 = pdf.get(5);
    const double f_6 = pdf.get(6);
    const double f_7 = pdf.get(7);
    const double f_8 = pdf.get(8);
    const double f_9 = pdf.get(9);
    const double f_10 = pdf.get(10);
    const double f_11 = pdf.get(11);
    const double f_12 = pdf.get(12);
    const double f_13 = pdf.get(13);
    const double f_14 = pdf.get(14);
    const double f_15 = pdf.get(15);
    const double f_16 = pdf.get(16);
    const double f_17 = pdf.get(17);
    const double f_18 = pdf.get(18);
    const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const double vel1Term = f_1 + f_11 + f_15 + f_7;
    const double vel2Term = f_12 + f_13 + f_5;
    const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
    out[offset] = rho;
  }
}

__global__ void kernel_set(
    cuda::FieldAccessor<double> pdf,
    const double *RESTRICT const rho_in) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, uint(1u));
    const double f_0 = pdf.get(0);
    const double f_1 = pdf.get(1);
    const double f_2 = pdf.get(2);
    const double f_3 = pdf.get(3);
    const double f_4 = pdf.get(4);
    const double f_5 = pdf.get(5);
    const double f_6 = pdf.get(6);
    const double f_7 = pdf.get(7);
    const double f_8 = pdf.get(8);
    const double f_9 = pdf.get(9);
    const double f_10 = pdf.get(10);
    const double f_11 = pdf.get(11);
    const double f_12 = pdf.get(12);
    const double f_13 = pdf.get(13);
    const double f_14 = pdf.get(14);
    const double f_15 = pdf.get(15);
    const double f_16 = pdf.get(16);
    const double f_17 = pdf.get(17);
    const double f_18 = pdf.get(18);
    const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
    const double momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
    const double vel1Term = f_1 + f_11 + f_15 + f_7;
    const double momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
    const double vel2Term = f_12 + f_13 + f_5;
    const double momdensity_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
    const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;

    // calculate current velocity (before density change)
    const double conversion = double(1) / rho;
    const double u_old[3] = {momdensity_0 * conversion, momdensity_1 * conversion, momdensity_2 * conversion};

    Equilibrium::kernel_set_device(pdf, u_old, rho_in[offset]);
  }
}

double get(
    cuda::GPUField<double> const *pdf_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(1u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  double rho = dev_data[0u];
  return rho;
}

void set(
    cuda::GPUField<double> *pdf_field,
    const double rho,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(1u, rho);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}

std::vector<double> get(
    cuda::GPUField<double> const *pdf_field,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(ci.numCells());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(dev_data_ptr);
  kernel();
  std::vector<double> out(ci.numCells());
  thrust::copy(dev_data.begin(), dev_data.end(), out.begin());
  return out;
}

void set(
    cuda::GPUField<double> *pdf_field,
    std::vector<double> const &values,
    CellInterval const &ci) {
  thrust::device_vector<double> dev_data(values.begin(), values.end());
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}
} // namespace Density

namespace Velocity {
__global__ void kernel_set(
    cuda::FieldAccessor<double> pdf,
    cuda::FieldAccessor<double> force,
    const double *RESTRICT const u_in) {
  pdf.set(blockIdx, threadIdx);
  force.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, uint(3u));
    const uint_t bufsize = 3u;
    const double *RESTRICT const u = u_in + bufsize * offset;
    const double f_0 = pdf.get(0);
    const double f_1 = pdf.get(1);
    const double f_2 = pdf.get(2);
    const double f_3 = pdf.get(3);
    const double f_4 = pdf.get(4);
    const double f_5 = pdf.get(5);
    const double f_6 = pdf.get(6);
    const double f_7 = pdf.get(7);
    const double f_8 = pdf.get(8);
    const double f_9 = pdf.get(9);
    const double f_10 = pdf.get(10);
    const double f_11 = pdf.get(11);
    const double f_12 = pdf.get(12);
    const double f_13 = pdf.get(13);
    const double f_14 = pdf.get(14);
    const double f_15 = pdf.get(15);
    const double f_16 = pdf.get(16);
    const double f_17 = pdf.get(17);
    const double f_18 = pdf.get(18);
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
    cuda::GPUField<double> *pdf_field,
    cuda::GPUField<double> *force_field,
    Vector3<double> const &u,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(u.data(), u.data() + 3u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_set);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*force_field, ci));
  kernel.addParam(const_cast<const double *>(dev_data_ptr));
  kernel();
}
} // namespace Velocity

namespace MomentumDensity {
__global__ void kernel_sum(
    cuda::FieldAccessor<double> pdf,
    cuda::FieldAccessor<double> force,
    double *RESTRICT const out) {
  pdf.set(blockIdx, threadIdx);
  force.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    const uint bufsize = 3u;
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, bufsize);
    const double f_0 = pdf.get(0);
    const double f_1 = pdf.get(1);
    const double f_2 = pdf.get(2);
    const double f_3 = pdf.get(3);
    const double f_4 = pdf.get(4);
    const double f_5 = pdf.get(5);
    const double f_6 = pdf.get(6);
    const double f_7 = pdf.get(7);
    const double f_8 = pdf.get(8);
    const double f_9 = pdf.get(9);
    const double f_10 = pdf.get(10);
    const double f_11 = pdf.get(11);
    const double f_12 = pdf.get(12);
    const double f_13 = pdf.get(13);
    const double f_14 = pdf.get(14);
    const double f_15 = pdf.get(15);
    const double f_16 = pdf.get(16);
    const double f_17 = pdf.get(17);
    const double f_18 = pdf.get(18);
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
    cuda::GPUField<double> const *pdf_field,
    cuda::GPUField<double> const *force_field) {
  thrust::device_vector<double> dev_data(3u, double{0});
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
    Cell cell(x, y, z);
    CellInterval ci(cell, cell);
    auto kernel = cuda::make_kernel(kernel_sum);
    kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
    kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*force_field, ci));
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
    cuda::FieldAccessor<double> pdf,
    double *RESTRICT const out) {
  pdf.set(blockIdx, threadIdx);
  if (pdf.isValidPosition()) {
    const uint bufsize = 9u;
    const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, bufsize);
    const double f_0 = pdf.get(0);
    const double f_1 = pdf.get(1);
    const double f_2 = pdf.get(2);
    const double f_3 = pdf.get(3);
    const double f_4 = pdf.get(4);
    const double f_5 = pdf.get(5);
    const double f_6 = pdf.get(6);
    const double f_7 = pdf.get(7);
    const double f_8 = pdf.get(8);
    const double f_9 = pdf.get(9);
    const double f_10 = pdf.get(10);
    const double f_11 = pdf.get(11);
    const double f_12 = pdf.get(12);
    const double f_13 = pdf.get(13);
    const double f_14 = pdf.get(14);
    const double f_15 = pdf.get(15);
    const double f_16 = pdf.get(16);
    const double f_17 = pdf.get(17);
    const double f_18 = pdf.get(18);
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
    cuda::GPUField<double> const *pdf_field,
    Cell const &cell) {
  CellInterval ci(cell, cell);
  thrust::device_vector<double> dev_data(9u);
  auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
  auto kernel = cuda::make_kernel(kernel_get);
  kernel.addFieldIndexingParam(cuda::FieldIndexing<double>::interval(*pdf_field, ci));
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