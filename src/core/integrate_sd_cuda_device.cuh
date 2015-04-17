#pragma once

typedef unsigned long long ull;
// atomicAdd implementation for double
__device__ double atomicAdd(double * address, double inc);

__device__ int fold_back_up(int value, int max);

__device__ int reduce_max(int * shared_cache, int value);

__device__ __inline__ real read_without_caching( const real * addr);

__device__ void reduce_sum(real * shared_cache);
