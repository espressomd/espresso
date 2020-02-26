#include "CudaHostAllocator.hpp"

void cuda_malloc_host(void **p, size_t n) {
  cudaError_t error = cudaMallocHost(p, n);

  if (error) {
    throw std::bad_alloc();
  }
}

void cuda_free_host(void *p) { cudaFreeHost(p); }
