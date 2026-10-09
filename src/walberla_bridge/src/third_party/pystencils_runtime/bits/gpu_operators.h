#pragma once

namespace pystencils::runtime::gpu
{
    struct Max
    {
        template <typename VT>
        __device__ __host__ auto operator()(const VT &a, const VT &b)
        {
            return a > b ? a : b;
        }
    };

    struct Min
    {
        template <typename VT>
        __device__ __host__ auto operator()(const VT &a, const VT &b)
        {
            return a < b ? a : b;
        }
    };

    struct Mul
    {
        template <typename VT>
        __device__ __host__ auto operator()(const VT &a, const VT &b)
        {
            return a * b;
        }
    };

    struct Div
    {
        template <typename VT>
        __device__ __host__ auto operator()(const VT &a, const VT &b)
        {
            return a / b;
        }
    };

    struct Add
    {
        template <typename VT>
        __device__ __host__ auto operator()(const VT &a, const VT &b)
        {
            return a + b;
        }
    };

} // pystencils::runtime::gpu