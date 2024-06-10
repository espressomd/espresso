/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "init.hpp"
#include "utils.cuh"

#include <utils/constants.hpp>

#include <cuda.h>
#include <cuda_runtime.h>

#include <cstring>
#include <string>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

#ifdef CUDA

/** \name minimally required compute capability. */
/**@{*/
static const int computeCapabilityMinMajor = 3;
static const int computeCapabilityMinMinor = 0;
/**@}*/

void cuda_init() { CUDA_CHECK(cudaStreamCreate(&stream[0])) }

int cuda_get_n_gpus() {
  int deviceCount;
  CUDA_CHECK(cudaGetDeviceCount(&deviceCount))
  return deviceCount;
}

int cuda_check_gpu_compute_capability(int dev) {
  cudaDeviceProp deviceProp;
  CUDA_CHECK(cudaGetDeviceProperties(&deviceProp, dev))
  if (deviceProp.major < computeCapabilityMinMajor ||
      (deviceProp.major == computeCapabilityMinMajor &&
       deviceProp.minor < computeCapabilityMinMinor)) {
    return ES_ERROR;
  }
  return ES_OK;
}

/**
 * @brief Safely copy the device name and pad the string with null characters.
 */
static void cuda_copy_gpu_name(char *const name, cudaDeviceProp const &prop) {
  char buffer[256] = {'\0'};
  std::strncpy(buffer, prop.name, 256);
  name[255] = '\0';
  std::strncpy(name, buffer, 256);
}

void cuda_get_gpu_name(int dev, char *const name) {
  cudaDeviceProp deviceProp;
  CUDA_CHECK(cudaGetDeviceProperties(&deviceProp, dev))
  cuda_copy_gpu_name(name, deviceProp);
}

EspressoGpuDevice cuda_get_device_props(const int dev) {
  cudaDeviceProp deviceProp;
  CUDA_CHECK(cudaGetDeviceProperties(&deviceProp, dev))
  EspressoGpuDevice device{dev,
                           "",
                           "",
                           -1,
                           deviceProp.major,
                           deviceProp.minor,
                           deviceProp.totalGlobalMem,
                           deviceProp.multiProcessorCount};
  cuda_copy_gpu_name(device.name, deviceProp);
  return device;
}

void cuda_set_device(int dev) {
  CUDA_CHECK(cudaSetDevice(dev))
  CUDA_CHECK(cudaStreamDestroy(stream[0]))
  CUDA_CHECK(cudaStreamCreate(&stream[0]))
}

int cuda_get_device() {
  int dev;
  CUDA_CHECK(cudaGetDevice(&dev))
  return dev;
}

int cuda_test_device_access() {
  int *d = nullptr;
  int h = 42;
  cudaError_t err;

  err = cudaMalloc((void **)&d, sizeof(int));
  if (err != cudaSuccess) {
    throw cuda_runtime_error_cuda(err);
  }
  err = cudaMemcpy(d, &h, sizeof(int), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cudaFree(d);
    throw cuda_runtime_error_cuda(err);
  }
  h = 0;
  err = cudaMemcpy(&h, d, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d);
  if (err != cudaSuccess) {
    throw cuda_runtime_error_cuda(err);
  }
  if (h != 42) {
    return ES_ERROR;
  }
  return ES_OK;
}

void cuda_check_device() {
  if (cuda_get_n_gpus() == 0) {
    throw cuda_runtime_error("No GPU was found.");
  }
  auto const devID = cuda_get_device();
  auto const compute_capability = cuda_check_gpu_compute_capability(devID);
  auto const communication_test = cuda_test_device_access();
  if (compute_capability != ES_OK or communication_test != ES_OK) {
    throw cuda_runtime_error("CUDA device " + std::to_string(devID) +
                             " is not capable of running ESPResSo.");
  }
}

#endif /* defined(CUDA) */
