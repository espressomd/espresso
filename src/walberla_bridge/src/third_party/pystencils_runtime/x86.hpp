/*
Copyright 2026 Frederik Hennig

This file is part of pystencils.

pystencils is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

pystencils is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with pystencils.
If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include "./bits/platforms.h"

#if PYSTENCILS_PLATFORM_HAS_SSE || PYSTENCILS_PLATFORM_HAS_AVX || PYSTENCILS_PLATFORM_HAS_AVX512
#include <immintrin.h>
#endif

#if PYSTENCILS_PLATFORM_HAS_SSE
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <nmmintrin.h>
#endif

#include "./bits/simd_horizontal_helpers.h"
#include "./bits/philox_rand.h"

/**
 * @brief Runtime library for the x86 code generation platform
 */

namespace pystencils::runtime::x86
{

    /**
     * SHORT ARRAYS
     */

#if PYSTENCILS_PLATFORM_HAS_SSE
    struct __m128dx2
    {
        __m128d v[2];

        __m128d &operator[](size_t i)
        {
            return v[i];
        }

        const __m128d &operator[](size_t i) const
        {
            return v[i];
        }
    };

    struct __m128x4
    {
        __m128 v[4];

        __m128 &operator[](size_t i)
        {
            return v[i];
        }

        const __m128 &operator[](size_t i) const
        {
            return v[i];
        }
    };
#endif

#if PYSTENCILS_PLATFORM_HAS_AVX
    struct __m256dx2
    {
        __m256d v[2];

        __m256d &operator[](size_t i)
        {
            return v[i];
        }

        const __m256d &operator[](size_t i) const
        {
            return v[i];
        }
    };

    struct __m256x4
    {
        __m256 v[4];

        __m256 &operator[](size_t i)
        {
            return v[i];
        }

        const __m256 &operator[](size_t i) const
        {
            return v[i];
        }
    };
#endif

#if PYSTENCILS_PLATFORM_HAS_AVX512
    struct __m512dx2
    {
        __m512d v[2];

        __m512d &operator[](size_t i)
        {
            return v[i];
        }

        const __m512d &operator[](size_t i) const
        {
            return v[i];
        }
    };

    struct __m512x4
    {
        __m512 v[4];

        __m512 &operator[](size_t i)
        {
            return v[i];
        }

        const __m512 &operator[](size_t i) const
        {
            return v[i];
        }
    };
#endif

    /**
     * CUSTOM INTRINSICS
     */

    namespace intrin
    {

#if PYSTENCILS_PLATFORM_HAS_SSE
        inline __m128i _mm_pack_2x64_epi32(__m128i a, __m128i b)
        {
            __m128 result = _mm_shuffle_ps(_mm_castsi128_ps(a),
                                             _mm_castsi128_ps(b), _MM_SHUFFLE(2, 0, 2, 0));
            return _mm_castps_si128(result);
        }
#endif

#if PYSTENCILS_PLATFORM_HAS_AVX
        //  From https://stackoverflow.com/questions/69408063/how-to-convert-int-64-to-int-32-with-avx-but-without-avx-512
        inline __m256i _mm256_pack_2x64_epi32(__m256i a, __m256i b)
        {
            __m256 combined = _mm256_shuffle_ps(_mm256_castsi256_ps(a),
                                                _mm256_castsi256_ps(b), _MM_SHUFFLE(2, 0, 2, 0));
            __m256d ordered = _mm256_permute4x64_pd(_mm256_castps_pd(combined), _MM_SHUFFLE(3, 1, 2, 0));
            return _mm256_castpd_si256(ordered);
        }
#endif
    } // namespace intrin

    /**
     * RANDOM NUMBER GENERATION
     *
     * Characteristics of vector RNGs:
     *
     *  - `fptype` is the scalar floating-point type of the RNG results
     *  - `vtype = fptype x vfactor` is the SIMD result type (with `vfactor` lanes)
     *  - `idx_type` is the scalar counter type, which must have *the same bit-width* as `fptype`
     *  - `vidx_type = idx_type x vfactor` is the vectorized counter type
     *  - `k` is the multiplicity of the generated random numbers
     *
     * Signature prototype:
     *
     * (vidx_type... counters, idx_type... keys) -> array< vtype, k >
     *
     *  - Counters get passed in as vectors of indices
     *  - Keys remain scalar arguments
     */

    namespace random
    {
#if PYSTENCILS_PLATFORM_HAS_SSE

        template <typename Idx>
        __m128x4 philox_float32(__m128i ctr0, __m128i ctr1, __m128i ctr2, __m128i ctr3, Idx key0, Idx key1)
        {
            __m128x4 result;
            ::pystencils::runtime::random::detail::philox_float4(
                ctr0, ctr1, ctr2, ctr3,                    //
                (uint32_t)key0, (uint32_t)key1,            //
                result[0], result[1], result[2], result[3] //
            );
            return result;
        }

        template <typename Idx>
        __m128dx2 philox_float64(__m128i ctr0, __m128i ctr1, __m128i ctr2, __m128i ctr3, Idx key0, Idx key1)
        {
            using namespace ::pystencils::runtime::x86::intrin;

            //  convert 64-bit integer counters to 32-bit
            __m128i ctr0_i32{_mm_pack_2x64_epi32(ctr0, _mm_setzero_si128())};
            __m128i ctr1_i32{_mm_pack_2x64_epi32(ctr1, _mm_setzero_si128())};
            __m128i ctr2_i32{_mm_pack_2x64_epi32(ctr2, _mm_setzero_si128())};
            __m128i ctr3_i32{_mm_pack_2x64_epi32(ctr3, _mm_setzero_si128())};

            __m128dx2 result;
            __m128d ignore;
            ::pystencils::runtime::random::detail::philox_double2(
                ctr0_i32, ctr1_i32, ctr2_i32, ctr3_i32, //
                (uint32_t)key0, (uint32_t)key1,         //
                result[0], ignore, result[1], ignore    //
            );
            return result;
        }
#endif

#if PYSTENCILS_PLATFORM_HAS_AVX
        template <typename Idx>
        __m256x4 philox_float32(__m256i ctr0, __m256i ctr1, __m256i ctr2, __m256i ctr3, Idx key0, Idx key1)
        {
            __m256x4 result;
            ::pystencils::runtime::random::detail::philox_float4(
                ctr0, ctr1, ctr2, ctr3,                    //
                (uint32_t)key0, (uint32_t)key1,            //
                result[0], result[1], result[2], result[3] //
            );
            return result;
        }

        template <typename Idx>
        __m256dx2 philox_float64(__m256i ctr0, __m256i ctr1, __m256i ctr2, __m256i ctr3, Idx key0, Idx key1)
        {
            using namespace ::pystencils::runtime::x86::intrin;

            //  convert 64-bit integer counters to 32-bit
            __m256i ctr0_i32{_mm256_pack_2x64_epi32(ctr0, _mm256_setzero_si256())};
            __m256i ctr1_i32{_mm256_pack_2x64_epi32(ctr1, _mm256_setzero_si256())};
            __m256i ctr2_i32{_mm256_pack_2x64_epi32(ctr2, _mm256_setzero_si256())};
            __m256i ctr3_i32{_mm256_pack_2x64_epi32(ctr3, _mm256_setzero_si256())};

            __m256dx2 result;
            __m256d ignore;
            ::pystencils::runtime::random::detail::philox_double2(
                ctr0_i32, ctr1_i32, ctr2_i32, ctr3_i32, //
                (uint32_t)key0, (uint32_t)key1,         //
                result[0], ignore, result[1], ignore    //
            );
            return result;
        }
#endif

#if PYSTENCILS_PLATFORM_HAS_AVX512
        template <typename Idx>
        __m512x4 philox_float32(__m512i ctr0, __m512i ctr1, __m512i ctr2, __m512i ctr3, Idx key0, Idx key1)
        {
            __m512x4 result;
            ::pystencils::runtime::random::detail::philox_float4(
                ctr0, ctr1, ctr2, ctr3,                    //
                (uint32_t)key0, (uint32_t)key1,            //
                result[0], result[1], result[2], result[3] //
            );
            return result;
        }

        template <typename Idx>
        __m512dx2 philox_float64(__m512i ctr0, __m512i ctr1, __m512i ctr2, __m512i ctr3, Idx key0, Idx key1)
        {
            //  convert 64-bit integer counters to 32-bit
            __m512i ctr0_i32{_mm512_zextsi256_si512(_mm512_cvtepi64_epi32(ctr0))};
            __m512i ctr1_i32{_mm512_zextsi256_si512(_mm512_cvtepi64_epi32(ctr1))};
            __m512i ctr2_i32{_mm512_zextsi256_si512(_mm512_cvtepi64_epi32(ctr2))};
            __m512i ctr3_i32{_mm512_zextsi256_si512(_mm512_cvtepi64_epi32(ctr3))};

            __m512dx2 result;
            __m512d ignore;
            ::pystencils::runtime::random::detail::philox_double2(
                ctr0_i32, ctr1_i32, ctr2_i32, ctr3_i32, //
                (uint32_t)key0, (uint32_t)key1,         //
                result[0], ignore, result[1], ignore    //
            );
            return result;
        }
#endif
    } // namespace random

} // namespace pystencils::runtime::x86
