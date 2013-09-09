/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  
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


#include "utils.hpp"
#include "cuda_init.hpp"


#ifdef CUDA

#include <cuda.h>

/** \name minimally required compute capability. */
/*@{*/
static const int computeCapabilityMinMajor = 1;
static const int computeCapabilityMinMinor = 1;
/*@}*/

const char *cuda_error;

/// get the number of CUDA devices.
int cuda_get_n_gpus()
{
  int deviceCount;
  cudaError_t error = cudaGetDeviceCount(&deviceCount);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    return -1;
  }
  return deviceCount;
}

int cuda_check_gpu(int dev)
{
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

void cuda_get_gpu_name(int dev, char name[64])
{
  cudaDeviceProp deviceProp;
  cudaError_t error = cudaGetDeviceProperties(&deviceProp, dev);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    strcpy(name, "no GPU");
    return;
  }
  strncpy(name, deviceProp.name, 63);
  name[63] = 0;
}

int cuda_set_device(int dev)
{
  cudaError_t error = cudaSetDevice(dev);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    return ES_ERROR;
  }
  else
    return ES_OK;
}

int cuda_get_device()
{
  int dev;
  cudaError_t error = cudaGetDevice(&dev);
  if (error != cudaSuccess) {
    cuda_error = cudaGetErrorString(error);
    return -1;
  }
  else
    return dev;
}

#endif /* defined(CUDA) */
