#ifndef CUDA_WRAPPER_HPP
#define CUDA_WRAPPER_HPP

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
