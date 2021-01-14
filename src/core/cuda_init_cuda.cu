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

#include <cuda.h>

#include "cuda_init.hpp"
#include "cuda_utils.cuh"

#include <utils/constants.hpp>

#include <cstring>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

#ifdef CUDA

/** \name minimally required compute capability. */
/**@{*/
static const int computeCapabilityMinMajor = 3;
static const int computeCapabilityMinMinor = 0;
/**@}*/

void cuda_init() {
  cudaError_t error = cudaStreamCreate(&stream[0]);
  if (error != cudaSuccess) {
    throw cuda_runtime_error_impl(error);
  }
}

/// get the number of CUDA devices.
int cuda_get_n_gpus() {
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);
  if (error != cudaSuccess) {
    throw cuda_runtime_error_impl(error);
  }
  return deviceCount;
}

int cuda_check_gpu_compute_capability(int dev) {
  cudaDeviceProp deviceProp;
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, dev);
  if (error != cudaSuccess) {
    throw cuda_runtime_error_impl(error);
  }
  if (deviceProp.major < computeCapabilityMinMajor ||
      (deviceProp.major == computeCapabilityMinMajor &&
       deviceProp.minor < computeCapabilityMinMinor)) {
    return ES_ERROR;
  }
  return ES_OK;
}

void cuda_get_gpu_name(int dev, char name[64]) {
  cudaDeviceProp deviceProp;
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, dev);
  if (error != cudaSuccess) {
    std::strncpy(name, "no GPU", 63);
    throw cuda_runtime_error_impl(error);
  }
  std::strncpy(name, deviceProp.name, 63);
  name[63] = 0;
}

EspressoGpuDevice cuda_get_device_props(const int dev) {
  cudaDeviceProp deviceProp;
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, dev);
  if (error != cudaSuccess) {
    throw cuda_runtime_error_impl(error);
  }
  EspressoGpuDevice device{dev,
                           "",
                           "",
                           -1,
                           deviceProp.major,
                           deviceProp.minor,
                           deviceProp.totalGlobalMem,
                           deviceProp.multiProcessorCount};
  std::strncpy(device.name, deviceProp.name, 64);
  device.name[63] = '\0';
  return device;
}

void cuda_set_device(int dev) {
  cudaError_t error = cudaSetDevice(dev);
  if (error != cudaSuccess) {
    throw cuda_runtime_error_impl(error);
  }
  error = cudaStreamDestroy(stream[0]);
  if (error != cudaSuccess) {
    throw cuda_runtime_error_impl(error);
  }
  error = cudaStreamCreate(&stream[0]);
  if (error != cudaSuccess) {
    throw cuda_runtime_error_impl(error);
  }
}

int cuda_get_device() {
  int dev;
  cudaError_t error = cudaGetDevice(&dev);
  if (error != cudaSuccess) {
    throw cuda_runtime_error_impl(error);
  }
  return dev;
}

int cuda_test_device_access() {
  int *d = nullptr;
  int h = 42;
  cudaError_t err;

  err = cudaMalloc((void **)&d, sizeof(int));
  if (err != cudaSuccess) {
    throw cuda_runtime_error_impl(err);
  }
  err = cudaMemcpy(d, &h, sizeof(int), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cudaFree(d);
    throw cuda_runtime_error_impl(err);
  }
  h = 0;
  err = cudaMemcpy(&h, d, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d);
  if (err != cudaSuccess) {
    throw cuda_runtime_error_impl(err);
  }
  if (h != 42) {
    return ES_ERROR;
  }
  return ES_OK;
}

#endif /* defined(CUDA) */
