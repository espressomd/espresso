#pragma once

namespace pystencils::runtime::gpu {

// No direct implementation for all atomic operations available
// -> add support by custom implementations using a CAS mechanism

// - atomicMul (double/float)
//   see https://stackoverflow.com/questions/43354798/atomic-multiplication-and-division
__device__ __forceinline__ double atomicMul(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int oldValue = *address_as_ull, assumed;
    do {
      assumed = oldValue;
      oldValue = atomicCAS(address_as_ull, assumed, __double_as_longlong(val *
                           __longlong_as_double(assumed)));
    } while (assumed != oldValue);

    return __longlong_as_double(oldValue);
}

__device__ __forceinline__ float atomicMul(float* address, float val) {
    int* address_as_int = (int*)address;
    int old = *address_as_int;
    int assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_int, assumed, __float_as_int(val * __int_as_float(assumed)));
    } while (assumed != old);

    return __int_as_float(old);
}

// - atomicMin (double/float)
//   see https://stackoverflow.com/questions/17399119/how-do-i-use-atomicmax-on-floating-point-values-in-cuda
__device__ __forceinline__ double atomicMin(double *address, double val)
{
    unsigned long long ret = __double_as_longlong(*address);
    while(val < __longlong_as_double(ret))
    {
        unsigned long long old = ret;
        if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val))) == old)
            break;
    }
    return __longlong_as_double(ret);
}

__device__ __forceinline__ float atomicMin(float *address, float val)
{
    int ret = __float_as_int(*address);
    while(val < __int_as_float(ret))
    {
        int old = ret;
        if((ret = atomicCAS((int *)address, old, __float_as_int(val))) == old)
            break;
    }
    return __int_as_float(ret);
}

// - atomicMax (double/float)
//   see https://stackoverflow.com/questions/17399119/how-do-i-use-atomicmax-on-floating-point-values-in-cuda
__device__ __forceinline__ double atomicMax(double *address, double val)
{
    unsigned long long ret = __double_as_longlong(*address);
    while(val > __longlong_as_double(ret))
    {
        unsigned long long old = ret;
        if((ret = atomicCAS((unsigned long long *)address, old, __double_as_longlong(val))) == old)
            break;
    }
    return __longlong_as_double(ret);
}

__device__ __forceinline__ float atomicMax(float *address, float val)
{
    int ret = __float_as_int(*address);
    while(val > __int_as_float(ret))
    {
        int old = ret;
        if((ret = atomicCAS((int *)address, old, __float_as_int(val))) == old)
            break;
    }
    return __int_as_float(ret);
}

} // pystencils::runtime::gpu
