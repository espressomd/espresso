// This works around a bug which prevents compilation of trhust with Clang
// https://github.com/thrust/thrust/issues/1032 (comment)

#ifndef CLANG_THRUST_WORKAROUND_CUH
#define CLANG_THRUST_WORKAROUND_CUH

#ifdef __clang__
#define THRUST_CUB_NS_PREFIX namespace thrust::cuda_cub {
#define THRUST_CUB_NS_POSTFIX }

#include <thrust/system/cuda/detail/cub/util_debug.cuh>

using namespace thrust::cuda_cub::cub;
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#endif

#endif
