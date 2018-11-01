#ifndef CURAND_WRAPPER_HPP
#define CURAND_WRAPPER_HPP

#if defined(__CUDACC__)

#include <curand_kernel.h>

#elif defined(__HIPCC__)

#include <rocrand/rocrand_kernel.h>

class philox4x32_10_stateless : private rocrand_device::philox4x32_10_engine {
public:
  FQUALIFIERS
  philox4x32_10_stateless() {}

  FQUALIFIERS
  uint4 operator()(uint4 counter, uint2 key) {
    return ten_rounds(counter, key);
  }
};

__forceinline__ __device__ uint4 curand_Philox4x32_10(uint4 counter,
                                                      uint2 key) {
  philox4x32_10_stateless *e;
  return (*e)(counter, key);
}

#endif

#endif // CURAND_WRAPPER_HPP
