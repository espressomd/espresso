#ifndef UTILS_DEVICE_QUALIFIER_HPP
#define UTILS_DEVICE_QUALIFIER_HPP

#if defined(__HIPCC__) or defined(__CUDACC__)
#define DEVICE_CODE
#endif

#ifdef DEVICE_CODE
#define DEVICE_QUALIFIER __host__ __device__
#define DEVICE_ASSERT(A) void((A))
#else
#define DEVICE_QUALIFIER
#define DEVICE_ASSERT(A) assert((A))
#endif

#endif // ESPRESSO_DEVICE_QUALIFIER_HPP
