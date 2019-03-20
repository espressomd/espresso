#ifndef UTILS_DEVICE_QUALIFIER_HPP
#define UTILS_DEVICE_QUALIFIER_HPP

#if defined(__CUDACC__)
#define DEVICE_QUALIFIER __host__ __device__
#else
#define DEVICE_QUALIFIER
#endif
#endif //ESPRESSO_DEVICE_QUALIFIER_HPP
