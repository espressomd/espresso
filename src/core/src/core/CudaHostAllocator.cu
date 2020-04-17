#include "CudaHostAllocator.hpp"

#include "cuda_wrapper.hpp"

#include <stdexcept>

void cuda_malloc_host(void **p, std::size_t n) {
  cudaError_t error = cudaMallocHost(p, n);

  if (error) {
    throw std::bad_alloc();
  }
}

void cuda_free_host(void *p) { cudaFreeHost(p); }
