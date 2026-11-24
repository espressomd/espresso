/*
 * Copyright (C) 2026 The ESPResSo project
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

#pragma once

/*
 * CUDA-to-HIP compatibility shim. When the GPU sources are compiled as HIP
 * (clang -x hip, selected by the ESPRESSO_BUILD_WITH_HIP CMake option) this
 * header includes the HIP runtime and hipFFT and aliases the CUDA spellings
 * used by the ESPResSo core GPU translation units to their HIP equivalents.
 * On NVIDIA it includes the CUDA runtime and cuFFT and is otherwise a no-op,
 * so the NVIDIA build is unchanged. This is the only file that maps CUDA names
 * to HIP.
 */

#if defined(__HIP__)

#include <cstdlib>
#include <cstring>

#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>

using cudaError_t = hipError_t;
using cudaStream_t = hipStream_t;
using cudaDeviceProp = hipDeviceProp_t;

#define cudaSuccess hipSuccess
#define cudaErrorInvalidValue hipErrorInvalidValue
#define cudaErrorInvalidDevice hipErrorInvalidDevice
#define cudaErrorNotPermitted hipErrorNotSupported
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost

#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError
#define cudaGetDevice hipGetDevice
#define cudaSetDevice hipSetDevice
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaStreamCreate hipStreamCreate
#define cudaStreamDestroy hipStreamDestroy
#define cudaMalloc hipMalloc
#define cudaFree hipFree
#define cudaMemcpy hipMemcpy
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemset hipMemset
#define cudaMemsetAsync hipMemsetAsync
#define cudaFreeHost hipHostFree

// hipHostMalloc is a template on the pointer type and takes a third flags
// argument, so it does not bind the two-argument void** call used by the
// pinned-host allocator (and the nullptr probe in the unit test) the way the
// non-template cudaMallocHost(void**, size_t) does. Forward explicitly.
static inline hipError_t cudaMallocHost(void **ptr, size_t size) {
  return hipHostMalloc(ptr, size, 0u);
}

using cufftComplex = hipfftComplex;
using cufftDoubleComplex = hipfftDoubleComplex;
using cufftHandle = hipfftHandle;

#define cufftExecR2C hipfftExecR2C
#define cufftExecC2R hipfftExecC2R
#define cufftExecD2Z hipfftExecD2Z
#define cufftExecZ2D hipfftExecZ2D
#define cufftPlan3d hipfftPlan3d
#define cufftDestroy hipfftDestroy
#define CUFFT_R2C HIPFFT_R2C
#define CUFFT_C2R HIPFFT_C2R
#define CUFFT_D2Z HIPFFT_D2Z
#define CUFFT_Z2D HIPFFT_Z2D
#define CUFFT_SUCCESS HIPFFT_SUCCESS

#else /* NVIDIA */

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>

#endif
