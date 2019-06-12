#ifndef UTILS_DEVICE_QUALIFIER_HPP
#define UTILS_DEVICE_QUALIFIER_HPP

#if defined(__CUDACC__)
#define DEVICE_COMPILER cudacc
#elif defined(__HIPCC__)
#define DEVICE_COMPILER hipcc
#endif

#if defined(DEVICE_COMPILER)
#define DEVICE_THROW(E)
#define DEVICE_QUALIFIER __host__ __device__
#define DEVICE_ASSERT(A) void((A))
#else
#define DEVICE_THROW(E) throw(E)
#define DEVICE_QUALIFIER
#define DEVICE_ASSERT(A) assert((A))
#endif

#endif // ESPRESSO_DEVICE_QUALIFIER_HPP
