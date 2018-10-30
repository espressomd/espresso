#ifndef CURAND_WRAPPER_HPP
#define CURAND_WRAPPER_HPP

#if defined(__HIPCC__) // AMD or Nvidia via HIP

#include <hiprand_kernel.h>

#else

#include <curand_kernel.h>

#endif
#endif
