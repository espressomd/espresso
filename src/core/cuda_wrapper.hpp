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
#ifndef CUDA_WRAPPER_HPP
#define CUDA_WRAPPER_HPP

#if defined(__CUDACC__)

#include <cuda.h>

#endif

#if defined(__HIPCC__) // AMD or Nvidia via HIP

#include <hip/hip_runtime.h>

#define cudaDeviceProp hipDeviceProp_t
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaErrorInvalidValue hipErrorInvalidValue
#define cudaError_t hipError_t
#define cudaEventCreate hipEventCreate
#define cudaEventDestroy hipEventDestroy
#define cudaEventElapsedTime hipEventElapsedTime
#define cudaEventRecord hipEventRecord
#define cudaEventSynchronize hipEventSynchronize
#define cudaEvent_t hipEvent_t
#define cudaFree hipFree
#define cudaFuncSetCacheConfig(a, b)
#define cudaGetDevice hipGetDevice
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError
#define cudaGetSymbolAddress hipGetSymbolAddress
#define cudaFreeHost hipHostFree
#define cudaHostAlloc hipHostMalloc
#define cudaHostAllocWriteCombined hipHostMallocWriteCombined
#define cudaMalloc hipMalloc
#define cudaMallocHost hipHostMalloc
#define cudaMemcpy hipMemcpy
#define cudaMemcpy2D hipMemcpy2D
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyToSymbol hipMemcpyToSymbol
#define cudaMemGetInfo hipMemGetInfo
#define cudaMemset hipMemset
#define cudaMemsetAsync hipMemsetAsync
#define cudaSetDevice hipSetDevice
#define cudaStreamCreate hipStreamCreate
#define cudaStreamDestroy hipStreamDestroy
#define cudaStream_t hipStream_t
#define cudaSuccess hipSuccess
#define __all_sync __all

#else // Nvidia via CUDA

#define hipLaunchKernelGGL(kernel, blocks, threads, mem, stream, ...)          \
  do {                                                                         \
    kernel<<<blocks, threads, mem, stream>>>(__VA_ARGS__);                     \
  } while (0)
#define HIP_DYNAMIC_SHARED(type, var) extern __shared__ type var[];
#define HIP_SYMBOL(X) X

#endif

#if defined(__HIPCC__) and not defined(__CUDACC__) // AMD-only

#define make_uint3 dim3

#endif

#endif // CUDA_WRAPPER_HPP
