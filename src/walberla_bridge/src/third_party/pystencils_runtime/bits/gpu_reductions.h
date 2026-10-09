#pragma once

#if defined(__HIPCC__)
    #include <hipcub/hipcub.hpp>

    namespace backend
    {
        using namespace hipcub;
    }
#elif defined(__CUDACC__)
    #include <cub/block/block_reduce.cuh>

    namespace backend
    {
        using namespace cub;
    }
#else
    #error "Unsupported compiler. Needs either CUDA or HIP"
#endif

namespace pystencils::runtime::gpu::cub
{
    template <typename VT, typename Op, int BX, int BY = 1, int BZ = 1>
    __device__ __forceinline__ VT block_reduce(VT value)
    {
        using BlockReduce = backend::BlockReduce<VT, BX, backend::BLOCK_REDUCE_RAKING_COMMUTATIVE_ONLY, BY, BZ>;

        __shared__ typename BlockReduce::TempStorage temp_storage;

        return BlockReduce(temp_storage).Reduce(value, Op{});
    }
}