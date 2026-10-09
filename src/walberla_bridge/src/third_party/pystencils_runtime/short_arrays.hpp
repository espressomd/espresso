#pragma once

#include "./bits/platforms.h"

#if !PYSTENCILS_PLATFORM_IS_GPU

#include <cstddef>

using std::size_t;

#endif

#if PYSTENCILS_PLATFORM_HAS_SSE || PYSTENCILS_PLATFORM_HAS_AVX || PYSTENCILS_PLATFORM_HAS_AVX512
#include <immintrin.h>
#endif


namespace pystencils::runtime
{

    template <typename T, size_t N>
    struct ShortArray
    {
        T v[4];

        PYSTENCILS_DEVICE_QUALIFIERS T &operator[](size_t i)
        {
            return v[i];
        }

        PYSTENCILS_DEVICE_QUALIFIERS const T &operator[](size_t i) const
        {
            return v[i];
        }
    };


} // namespace pystencils::runtime

#undef QUALIFIERS
