#ifndef CORE_LB_GPU_EXPORT_HPP
#define CORE_LB_GPU_EXPORT_HPP

#ifdef __CUDACC__
#define ESPRESSO_GPU_EXPORT __device__
#else
#define ESPRESSO_GPU_EXPORT
#endif
#endif
