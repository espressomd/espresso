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

#include "cuda_init.hpp"
#include "cuda_utils.hpp"
#include "debug.hpp"

#include <utils/constants.hpp>

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

#ifdef CUDA

/** \name minimally required compute capability. */
/*@{*/
static const int computeCapabilityMinMajor = 3;
static const int computeCapabilityMinMinor = 0;
/*@}*/

const char *cuda_error;

void cuda_init() {
#if defined(__HIPCC__) and not defined(__CUDACC__)
  setenv("HSA_ENABLE_INTERRUPT", "0", 1);
#endif
#if defined(__HIPCC__) and not defined(__CUDACC__) and                         \
    HIP_VERSION_PATCH <= 19171 /* i.e. <= v2.4.0 */
  // Catch an exception that causes `import espressomd` to crash in
  // Python when no compatible GPU is available (fixed in HIP v2.5.0)
  try {
    cudaStreamCreate(&stream[0]);
  } catch (std::exception &) {
    // The ihipException is not part of the public HIP library,
    // so we catch its parent class std::exception instead
  }
#else
  cudaStreamCreate(&stream[0]);
#endif
}

/// get the number of CUDA devices.
int cuda_get_n_gpus() {
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    return -1;
  }
  return deviceCount;
}

int cuda_check_gpu(int dev) {
  cudaDeviceProp deviceProp;
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, dev);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    return ES_ERROR;
  }
  if (deviceProp.major < computeCapabilityMinMajor ||
      (deviceProp.major == computeCapabilityMinMajor &&
       deviceProp.minor < computeCapabilityMinMinor)) {
    cuda_error = "compute capability insufficient";
    return ES_ERROR;
  }
  return ES_OK;
}

void cuda_get_gpu_name(int dev, char name[64]) {
  cudaDeviceProp deviceProp;
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, dev);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    strncpy(name, "no GPU", 63);
  } else {
    strncpy(name, deviceProp.name, 63);
  }
  name[63] = 0;
}

int cuda_get_device_props(const int dev, EspressoGpuDevice &d) {
  cudaDeviceProp deviceProp;
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, dev);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    return ES_ERROR;
  }
  strncpy(d.name, deviceProp.name, 64);
  d.id = dev;
  d.total_memory = deviceProp.totalGlobalMem;
  d.compute_capability_major = deviceProp.major;
  d.compute_capability_minor = deviceProp.minor;
  d.n_cores = deviceProp.multiProcessorCount;

  return ES_OK;
}

int cuda_set_device(int dev) {
  cudaSetDevice(dev);
  cudaStreamDestroy(stream[0]);
  cudaError_t error = cudaStreamCreate(&stream[0]);

  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    throw std::runtime_error(cuda_error);
  }

  return ES_OK;
}

int cuda_get_device() {
  int dev;
  cudaError_t error = cudaGetDevice(&dev);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    return -1;
  }
  return dev;
}

int cuda_test_device_access() {
  int *d = nullptr;
  int h = 42;
  cudaError_t err;

  err = cudaMalloc((void **)&d, sizeof(int));
  if (err != cudaSuccess) {
    cuda_error = cudaGetErrorString(err);
    return ES_ERROR;
  }
  err = cudaMemcpy(d, &h, sizeof(int), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    cuda_error = cudaGetErrorString(err);
    return ES_ERROR;
  }
  h = 0;
  err = cudaMemcpy(&h, d, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d);

  if ((h == 42) && (err == cudaSuccess)) {
    return ES_OK;
  }
  cuda_error = cudaGetErrorString(err);
  return ES_ERROR;
}

#endif /* defined(CUDA) */
