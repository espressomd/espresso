/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "cuda_wrapper.hpp"

// TODO: throw exceptions upon errors initialization

#include "cuda_utils.hpp"
#include "cufft_wrapper.hpp"
#include "grid_based_algorithms/fd-electrostatics.cuh"
#include <stdexcept>
#include <string>
//#include <cuda_interface.hpp>
#include <cstdio>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

__global__ void createGreensfcn();
__global__ void multiplyGreensfcn(cufftComplex *charge_potential);

__device__ __constant__ FdElectrostatics::Parameters fde_parameters_gpu[1];

__device__ unsigned int fde_getThreadIndex() {

  return blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x +
         threadIdx.x;
}

__device__ cufftReal fde_getNode(int x, int y, int z) {
  auto *field =
      reinterpret_cast<cufftReal *>(fde_parameters_gpu->charge_potential);
  return field[fde_parameters_gpu->dim_y * fde_parameters_gpu->dim_x_padded *
                   z +
               fde_parameters_gpu->dim_x_padded * y + x];
}

__device__ void fde_setNode(int x, int y, int z, cufftReal value) {
  auto *field =
      reinterpret_cast<cufftReal *>(fde_parameters_gpu->charge_potential);
  field[fde_parameters_gpu->dim_y * fde_parameters_gpu->dim_x_padded * z +
        fde_parameters_gpu->dim_x_padded * y + x] = value;
}

__device__ cufftReal fde_getNode(int i) {
  int x = i % fde_parameters_gpu->dim_x_padded;
  i /= fde_parameters_gpu->dim_x_padded;
  int y = i % fde_parameters_gpu->dim_y;
  int z = i / fde_parameters_gpu->dim_y;
  return fde_getNode(x, y, z);
}

__device__ void fde_setNode(int i, cufftReal value) {
  int x = i % fde_parameters_gpu->dim_x_padded;
  i /= fde_parameters_gpu->dim_x_padded;
  int y = i % fde_parameters_gpu->dim_y;
  int z = i / fde_parameters_gpu->dim_y;
  fde_setNode(x, y, z, value);
}

FdElectrostatics::~FdElectrostatics() {
  cufftDestroy(plan_ifft);
  cufftDestroy(plan_fft);

  cuda_safe_mem(cudaFree(parameters.greensfcn));
  cuda_safe_mem(cudaFree(parameters.charge_potential));
}

FdElectrostatics::FdElectrostatics(InputParameters inputParameters,
                                   cudaStream_t stream)
    : parameters(inputParameters), cuda_stream(stream) {
  cuda_safe_mem(cudaMalloc((void **)&parameters.charge_potential,
                           sizeof(cufftComplex) * parameters.dim_z *
                               parameters.dim_y * (parameters.dim_x / 2 + 1)));

  cuda_safe_mem(cudaMalloc((void **)&parameters.greensfcn,
                           sizeof(cufftReal) * parameters.dim_z *
                               parameters.dim_y * (parameters.dim_x / 2 + 1)));

  if (cudaGetLastError() != cudaSuccess) {
    throw std::runtime_error("Failed to allocate");
  }

  cuda_safe_mem(cudaMemcpyToSymbol(HIP_SYMBOL(fde_parameters_gpu), &parameters,
                                   sizeof(Parameters)));

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (parameters.dim_z * parameters.dim_y * (parameters.dim_x / 2 + 1) +
       threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
  KERNELCALL_stream(createGreensfcn, dim_grid, threads_per_block, stream);

  /* create 3D FFT plans */

  if (cufftPlan3d(&plan_fft, parameters.dim_z, parameters.dim_y,
                  parameters.dim_x, CUFFT_R2C) != CUFFT_SUCCESS) {
    throw std::runtime_error("Unable to create fft plan");
  }

  if (cufftSetStream(plan_fft, cuda_stream) != CUFFT_SUCCESS) {
    throw std::runtime_error("Unable to assign FFT to cuda stream");
  }

  if (cufftPlan3d(&plan_ifft, parameters.dim_z, parameters.dim_y,
                  parameters.dim_x, CUFFT_C2R) != CUFFT_SUCCESS) {
    throw std::runtime_error("Unable to create ifft plan");
  }

  if (cufftSetStream(plan_ifft, cuda_stream) != CUFFT_SUCCESS) {
    throw std::runtime_error("Unable to assign FFT to cuda stream");
  }

  initialized = true;
}

__global__ void createGreensfcn() {
  unsigned int index = fde_getThreadIndex();
  unsigned int tmp;
  unsigned int coord[3];

  coord[0] = index % (fde_parameters_gpu->dim_x / 2 + 1);
  tmp = index / (fde_parameters_gpu->dim_x / 2 + 1);
  coord[1] = tmp % fde_parameters_gpu->dim_y;
  coord[2] = tmp / fde_parameters_gpu->dim_y;

  if (index < fde_parameters_gpu->dim_z * fde_parameters_gpu->dim_y *
                  (fde_parameters_gpu->dim_x / 2 + 1)) {

    if (index == 0) {
      // setting 0th Fourier mode to 0 enforces charge neutrality
      fde_parameters_gpu->greensfcn[index] = 0.0f;
    } else {
      fde_parameters_gpu->greensfcn[index] =
          -4.0f * PI_FLOAT * fde_parameters_gpu->prefactor *
          fde_parameters_gpu->agrid * fde_parameters_gpu->agrid * 0.5f /
          (cos(2.0f * PI_FLOAT * static_cast<cufftReal>(coord[0]) /
               static_cast<cufftReal>(fde_parameters_gpu->dim_x)) +
           cos(2.0f * PI_FLOAT * static_cast<cufftReal>(coord[1]) /
               static_cast<cufftReal>(fde_parameters_gpu->dim_y)) +
           cos(2.0f * PI_FLOAT * static_cast<cufftReal>(coord[2]) /
               static_cast<cufftReal>(fde_parameters_gpu->dim_z)) -
           3.0f) /
          static_cast<cufftReal>(fde_parameters_gpu->dim_x *
                                 fde_parameters_gpu->dim_y *
                                 fde_parameters_gpu->dim_z);
    }

    // fde_parameters_gpu->greensfcn[index] = 0.0f; //TODO delete
  }
}

__global__ void multiplyGreensfcn(cufftComplex *charge_potential) {

  unsigned int index = fde_getThreadIndex();

  if (index < fde_parameters_gpu->dim_z * fde_parameters_gpu->dim_y *
                  (fde_parameters_gpu->dim_x / 2 + 1)) {
    charge_potential[index].x *= fde_parameters_gpu->greensfcn[index];
    charge_potential[index].y *= fde_parameters_gpu->greensfcn[index];
  }
}

void FdElectrostatics::calculatePotential() {
  calculatePotential(parameters.charge_potential);
}

void FdElectrostatics::calculatePotential(cufftComplex *charge_potential) {

  if (cufftExecR2C(plan_fft, (cufftReal *)charge_potential, charge_potential) !=
      CUFFT_SUCCESS) {

    fprintf(stderr, "ERROR: Unable to execute FFT plan\n");
  }

  int threads_per_block = 64;
  int blocks_per_grid_y = 4;
  int blocks_per_grid_x =
      (parameters.dim_z * parameters.dim_y * (parameters.dim_x / 2 + 1) +
       threads_per_block * blocks_per_grid_y - 1) /
      (threads_per_block * blocks_per_grid_y);
  dim3 dim_grid = make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);

  KERNELCALL(multiplyGreensfcn, dim_grid, threads_per_block, charge_potential);

  if (cufftExecC2R(plan_ifft, charge_potential,
                   (cufftReal *)charge_potential) != CUFFT_SUCCESS) {

    fprintf(stderr, "ERROR: Unable to execute iFFT plan\n");
  }
}

FdElectrostatics::Grid FdElectrostatics::getGrid() {
  Grid g = {(float *)parameters.charge_potential, parameters.dim_x,
            parameters.dim_y, parameters.dim_z, parameters.agrid};
  return g;
}
