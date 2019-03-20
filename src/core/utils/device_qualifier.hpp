#ifndef UTILS_DEVICE_QUALIFIER_HPP
#define UTILS_DEVICE_QUALIFIER_HPP

#if defined(__CUDACC__)
#define DEVICE_QUALIFIER __host__ __device__
#else
#define DEVICE_QUALIFIER
#endif

#ifndef __CUDACC__
#define DEVICE_ASSERT(A) assert((A))
#else
#define DEVICE_ASSERT(A) void((A))
#endif

#endif // ESPRESSO_DEVICE_QUALIFIER_HPP
