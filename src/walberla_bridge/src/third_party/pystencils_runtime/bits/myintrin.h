/*
Copyright 2019-2023, Michael Kuron.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#if defined(__SSE2__) || (defined(_MSC_VER) && !defined(_M_ARM64))
QUALIFIERS __m128 _my_cvtepu32_ps(const __m128i v)
{
#if defined(__AVX512VL__) || defined(__AVX10_1__)
    return _mm_cvtepu32_ps(v);
#else
    __m128i v2 = _mm_srli_epi32(v, 1);
    __m128i v1 = _mm_and_si128(v, _mm_set1_epi32(1));
    __m128 v2f = _mm_cvtepi32_ps(v2);
    __m128 v1f = _mm_cvtepi32_ps(v1);
    return _mm_add_ps(_mm_add_ps(v2f, v2f), v1f);
#endif
}
#endif

#if defined(__SSE4_1__) || (defined(_MSC_VER) && !defined(_M_ARM64))
#if !defined(__AVX512VL__) && !defined(__AVX10_1__) && defined(__GNUC__) && __GNUC__ >= 5 && !defined(__clang__)
__attribute__((optimize("no-associative-math")))
#endif
QUALIFIERS __m128d _my_cvtepu64_pd(const __m128i x)
{
#if defined(__AVX512VL__) || defined(__AVX10_1__)
    return _mm_cvtepu64_pd(x);
#elif defined(__clang__)
    return __builtin_convertvector((uint64_t __attribute__((__vector_size__(16)))) x, __m128d);
#else
    __m128i xH = _mm_srli_epi64(x, 32);
    xH = _mm_or_si128(xH, _mm_castpd_si128(_mm_set1_pd(19342813113834066795298816.)));          //  2^84
    __m128i xL = _mm_blend_epi16(x, _mm_castpd_si128(_mm_set1_pd(0x0010000000000000)), 0xcc);   //  2^52
    __m128d f = _mm_sub_pd(_mm_castsi128_pd(xH), _mm_set1_pd(19342813118337666422669312.));     //  2^84 + 2^52
    return _mm_add_pd(f, _mm_castsi128_pd(xL));
#endif
}
#endif

#ifdef __AVX2__
QUALIFIERS __m256i _my256_set_m128i(__m128i hi, __m128i lo)
{
#if (!defined(__GNUC__) || __GNUC__ >= 8) || defined(__clang__)
    return _mm256_set_m128i(hi, lo);
#else
    return _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1);
#endif
}

QUALIFIERS __m256d _my256_set_m128d(__m128d hi, __m128d lo)
{
#if (!defined(__GNUC__) || __GNUC__ >= 8) || defined(__clang__)
    return _mm256_set_m128d(hi, lo);
#else
    return _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), hi, 1);
#endif
}

QUALIFIERS __m256 _my256_cvtepu32_ps(const __m256i v)
{
#if defined(__AVX512VL__) || defined(__AVX10_1__)
    return _mm256_cvtepu32_ps(v);
#else
    __m256i v2 = _mm256_srli_epi32(v, 1);
    __m256i v1 = _mm256_and_si256(v, _mm256_set1_epi32(1));
    __m256 v2f = _mm256_cvtepi32_ps(v2);
    __m256 v1f = _mm256_cvtepi32_ps(v1);
    return _mm256_add_ps(_mm256_add_ps(v2f, v2f), v1f);
#endif
}

#if !defined(__AVX512VL__) && !defined(__AVX10_1__) && defined(__GNUC__) && __GNUC__ >= 5 && !defined(__clang__)
__attribute__((optimize("no-associative-math")))
#endif
QUALIFIERS __m256d _my256_cvtepu64_pd(const __m256i x)
{
#if defined(__AVX512VL__) || defined(__AVX10_1__)
    return _mm256_cvtepu64_pd(x);
#elif defined(__clang__)
    return __builtin_convertvector((uint64_t __attribute__((__vector_size__(32)))) x, __m256d);
#else
    __m256i xH = _mm256_srli_epi64(x, 32);
    xH = _mm256_or_si256(xH, _mm256_castpd_si256(_mm256_set1_pd(19342813113834066795298816.)));          //  2^84
    __m256i xL = _mm256_blend_epi16(x, _mm256_castpd_si256(_mm256_set1_pd(0x0010000000000000)), 0xcc);   //  2^52
    __m256d f = _mm256_sub_pd(_mm256_castsi256_pd(xH), _mm256_set1_pd(19342813118337666422669312.));     //  2^84 + 2^52
    return _mm256_add_pd(f, _mm256_castsi256_pd(xL));
#endif
}
#endif

#if defined(__AVX512F__) || defined(__AVX10_512BIT__)
QUALIFIERS __m512i _my512_set_m128i(__m128i d, __m128i c, __m128i b, __m128i a)
{
    return _mm512_inserti32x4(_mm512_inserti32x4(_mm512_inserti32x4(_mm512_castsi128_si512(a), b, 1), c, 2), d, 3);
}

QUALIFIERS __m512d _my512_set_m256d(__m256d b, __m256d a)
{
    return _mm512_insertf64x4(_mm512_castpd256_pd512(a), b, 1);
}
#endif
