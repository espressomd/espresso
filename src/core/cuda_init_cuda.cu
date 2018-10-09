/*
  Copyright (C) 2010-2018 The ESPResSo project

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <hip/hip_runtime.h>

#include "cuda_init.hpp"
#include "cuda_utils.hpp"
#include "debug.hpp"
#include "utils.hpp"

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

#ifdef CUDA

/** \name minimally required compute capability. */
/*@{*/
static const int computeCapabilityMinMajor = 1;
static const int computeCapabilityMinMinor = 1;
/*@}*/

const char *cuda_error;

void cuda_init() { hipStreamCreate(&stream[0]); }

/// get the number of CUDA devices.
int cuda_get_n_gpus() {
  int deviceCount;
  hipError_t error = hipGetDeviceCount(&deviceCount);
  if (error != hipSuccess) {
    cuda_error = hipGetErrorString(error);
    return -1;
  }
  return deviceCount;
}

int cuda_check_gpu(int dev) {
  hipDeviceProp_t deviceProp;
  hipError_t error = hipGetDeviceProperties(&deviceProp, dev);
  if (error != hipSuccess) {
    cuda_error = hipGetErrorString(error);
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
  hipDeviceProp_t deviceProp;
  hipError_t error = hipGetDeviceProperties(&deviceProp, dev);
  if (error != hipSuccess) {
    cuda_error = hipGetErrorString(error);
    strcpy(name, "no GPU");
    return;
  }
  strncpy(name, deviceProp.name, 63);
  name[63] = 0;
}

int cuda_get_device_props(const int dev, EspressoGpuDevice &d) {
  hipDeviceProp_t deviceProp;
  hipError_t error = hipGetDeviceProperties(&deviceProp, dev);
  if (error != hipSuccess) {
    cuda_error = hipGetErrorString(error);
    return ES_ERROR;
  }
  strncpy(d.name, deviceProp.name, 64);
  d.id = dev;
  d.total_memory = deviceProp.totalGlobalMem;
  d.node = this_node;
  d.compute_capability_major = deviceProp.major;
  d.compute_capability_minor = deviceProp.minor;
  d.n_cores = deviceProp.multiProcessorCount;

  return ES_OK;
}

int cuda_set_device(int dev) {
  hipSetDevice(dev);
  hipStreamDestroy(stream[0]);
  hipError_t error = hipStreamCreate(&stream[0]);

  if (error != hipSuccess) {
    cuda_error = hipGetErrorString(error);
    throw std::runtime_error(cuda_error);
  }

  return ES_OK;
}

int cuda_get_device() {
  int dev;
  hipError_t error = hipGetDevice(&dev);
  if (error != hipSuccess) {
    cuda_error = hipGetErrorString(error);
    return -1;
  } else
    return dev;
}

int cuda_test_device_access() {
  int *d = 0;
  int h = 42;
  hipError_t err;

  err = hipMalloc((void **)&d, sizeof(int));
  if (err != hipSuccess) {
    cuda_error = hipGetErrorString(err);
    return ES_ERROR;
  }
  err = hipMemcpy(d, &h, sizeof(int), hipMemcpyHostToDevice);
  if (err != hipSuccess) {
    cuda_error = hipGetErrorString(err);
    return ES_ERROR;
  }
  h = 0;
  err = hipMemcpy(&h, d, sizeof(int), hipMemcpyDeviceToHost);
  hipFree(d);

  if ((h == 42) && (err == hipSuccess))
    return ES_OK;
  else {
    cuda_error = hipGetErrorString(err);
    return ES_ERROR;
  }
}

#endif /* defined(CUDA) */
