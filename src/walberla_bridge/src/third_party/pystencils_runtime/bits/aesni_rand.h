/*
Copyright 2010-2011, D. E. Shaw Research. All rights reserved.
Copyright 2019-2025, Michael Kuron.

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

#if defined(_MSC_VER) && defined(_M_ARM64)
#define __ARM_NEON
#endif

#ifdef __ARM_NEON
#include <arm_neon.h>
#if defined(__ARM_FEATURE_SVE)
#include <arm_sve.h>
#include <arm_neon_sve_bridge.h>
#endif
#else
#include <emmintrin.h> // SSE2
#include <wmmintrin.h> // AES
#ifdef __AVX__
#include <immintrin.h> // AVX*
#else
#include <smmintrin.h>  // SSE4
#ifdef __FMA__
#include <immintrin.h> // FMA
#endif
#endif
#endif
#include <cstdint>
#include <array>
#include <map>

#define QUALIFIERS inline
#define TWOPOW53_INV_DOUBLE (1.1102230246251565e-16)
#define TWOPOW32_INV_FLOAT (2.3283064e-10f)

#include "./myintrin.h"

typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

template <typename T, std::size_t Alignment>
class AlignedAllocator
{
public:
    typedef T value_type;

    template <typename U>
    struct rebind {
        typedef AlignedAllocator<U, Alignment> other;
    };

    T * allocate(const std::size_t n) const {
        if (n == 0) {
            return nullptr;
        }
#ifdef _WIN32
        void * const p = _aligned_malloc(n*sizeof(T), Alignment);
#else
        void * p;
        if (posix_memalign(&p, Alignment, n*sizeof(T)) != 0) {
          p = nullptr;
        }
#endif
        if (p == nullptr) {
            throw std::bad_alloc();
        }
        return static_cast<T *>(p);
    }

    void deallocate(T * const p, const std::size_t n) const {
#ifdef _WIN32
        _aligned_free(p);
#else
        free(p);
#endif
    }
};

template <typename Key, typename T>
using AlignedMap = std::map<Key, T, std::less<Key>, AlignedAllocator<std::pair<const Key, T>, sizeof(Key)>>;

#if defined(__AES__) || defined(_MSC_VER)
QUALIFIERS __m128i aesni_keygen_assist(__m128i temp1, __m128i temp2) {
    __m128i temp3; 
    temp2 = _mm_shuffle_epi32(temp2, 0xff); 
    temp3 = _mm_slli_si128(temp1, 0x4);
    temp1 = _mm_xor_si128(temp1, temp3);
    temp3 = _mm_slli_si128(temp3, 0x4);
    temp1 = _mm_xor_si128(temp1, temp3);
    temp3 = _mm_slli_si128(temp3, 0x4);
    temp1 = _mm_xor_si128(temp1, temp3);
    temp1 = _mm_xor_si128(temp1, temp2); 
    return temp1; 
}

QUALIFIERS std::array<__m128i,11> aesni_keygen(__m128i k) {
    std::array<__m128i,11> rk;
    __m128i tmp;
    
    rk[0] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x1);
    k = aesni_keygen_assist(k, tmp);
    rk[1] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x2);
    k = aesni_keygen_assist(k, tmp);
    rk[2] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x4);
    k = aesni_keygen_assist(k, tmp);
    rk[3] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x8);
    k = aesni_keygen_assist(k, tmp);
    rk[4] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x10);
    k = aesni_keygen_assist(k, tmp);
    rk[5] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x20);
    k = aesni_keygen_assist(k, tmp);
    rk[6] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x40);
    k = aesni_keygen_assist(k, tmp);
    rk[7] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x80);
    k = aesni_keygen_assist(k, tmp);
    rk[8] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x1b);
    k = aesni_keygen_assist(k, tmp);
    rk[9] = k;
    
    tmp = _mm_aeskeygenassist_si128(k, 0x36);
    k = aesni_keygen_assist(k, tmp);
    rk[10] = k;
    
    return rk;
}

QUALIFIERS const std::array<__m128i,11> & aesni_roundkeys(const __m128i & k128) {
    alignas(16) std::array<uint32,4> a;
    _mm_store_si128((__m128i*) a.data(), k128);
    
    static AlignedMap<std::array<uint32,4>, std::array<__m128i,11>> roundkeys;
    
    if(roundkeys.find(a) == roundkeys.end()) {
        auto rk = aesni_keygen(k128);
        roundkeys[a] = rk;
    }
    return roundkeys[a];
}

QUALIFIERS __m128i aesni1xm128i(const __m128i & in, const __m128i & k0) {
    auto k = aesni_roundkeys(k0);
    __m128i x = _mm_xor_si128(k[0], in);
    x = _mm_aesenc_si128(x, k[1]);
    x = _mm_aesenc_si128(x, k[2]);
    x = _mm_aesenc_si128(x, k[3]);
    x = _mm_aesenc_si128(x, k[4]);
    x = _mm_aesenc_si128(x, k[5]);
    x = _mm_aesenc_si128(x, k[6]);
    x = _mm_aesenc_si128(x, k[7]);
    x = _mm_aesenc_si128(x, k[8]);
    x = _mm_aesenc_si128(x, k[9]);
    x = _mm_aesenclast_si128(x, k[10]);
    return x;
}


QUALIFIERS void aesni_double2(uint32 ctr0, uint32 ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              double & rnd1, double & rnd2)
{
    // pack input and call AES
    __m128i c128 = _mm_set_epi32(ctr3, ctr2, ctr1, ctr0);
    __m128i k128 = _mm_set_epi32(key3, key2, key1, key0);
    c128 = aesni1xm128i(c128, k128);

    // convert 32 to 64 bit and put 0th and 2nd element into x, 1st and 3rd element into y
    __m128i x = _mm_and_si128(c128, _mm_set_epi32(0, 0xffffffff, 0, 0xffffffff));
    __m128i y = _mm_and_si128(c128, _mm_set_epi32(0xffffffff, 0, 0xffffffff, 0));
    y = _mm_srli_si128(y, 4);

    // calculate z = x ^ y << (53 - 32))
    __m128i z = _mm_sll_epi64(y, _mm_set1_epi64x(53 - 32));
    z = _mm_xor_si128(x, z);

    // convert uint64 to double
    __m128d rs = _my_cvtepu64_pd(z);
    // calculate rs * TWOPOW53_INV_DOUBLE + (TWOPOW53_INV_DOUBLE/2.0)
#ifdef __FMA__
    rs = _mm_fmadd_pd(rs, _mm_set1_pd(TWOPOW53_INV_DOUBLE), _mm_set1_pd(TWOPOW53_INV_DOUBLE/2.0));
#else
    rs = _mm_mul_pd(rs, _mm_set1_pd(TWOPOW53_INV_DOUBLE));
    rs = _mm_add_pd(rs, _mm_set1_pd(TWOPOW53_INV_DOUBLE/2.0));
#endif

    // store result
    alignas(16) double rr[2];
    _mm_store_pd(rr, rs);
    rnd1 = rr[0];
    rnd2 = rr[1];
}


QUALIFIERS void aesni_float4(uint32 ctr0, uint32 ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             float & rnd1, float & rnd2, float & rnd3, float & rnd4)
{
    // pack input and call AES
    __m128i c128 = _mm_set_epi32(ctr3, ctr2, ctr1, ctr0);
    __m128i k128 = _mm_set_epi32(key3, key2, key1, key0);
    c128 = aesni1xm128i(c128, k128);

    // convert uint32 to float
    __m128 rs = _my_cvtepu32_ps(c128);
    // calculate rs * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f)
#ifdef __FMA__
    rs = _mm_fmadd_ps(rs, _mm_set1_ps(TWOPOW32_INV_FLOAT), _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
#else
    rs = _mm_mul_ps(rs, _mm_set1_ps(TWOPOW32_INV_FLOAT));
    rs = _mm_add_ps(rs, _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
#endif

    // store result
    alignas(16) float r[4];
    _mm_store_ps(r, rs);
    rnd1 = r[0];
    rnd2 = r[1];
    rnd3 = r[2];
    rnd4 = r[3];
}


template<bool high>
QUALIFIERS __m128d _uniform_double_hq(__m128i x, __m128i y)
{
    // convert 32 to 64 bit
    if (high)
    {
        x = _mm_unpackhi_epi32(x, _mm_setzero_si128());
        y = _mm_unpackhi_epi32(y, _mm_setzero_si128());
    }
    else
    {
        x = _mm_unpacklo_epi32(x, _mm_setzero_si128());
        y = _mm_unpacklo_epi32(y, _mm_setzero_si128());
    }

    // calculate z = x ^ y << (53 - 32))
    __m128i z = _mm_sll_epi64(y, _mm_set1_epi64x(53 - 32));
    z = _mm_xor_si128(x, z);

    // convert uint64 to double
    __m128d rs = _my_cvtepu64_pd(z);
    // calculate rs * TWOPOW53_INV_DOUBLE + (TWOPOW53_INV_DOUBLE/2.0)
#ifdef __FMA__
    rs = _mm_fmadd_pd(rs, _mm_set1_pd(TWOPOW53_INV_DOUBLE), _mm_set1_pd(TWOPOW53_INV_DOUBLE/2.0));
#else
    rs = _mm_mul_pd(rs, _mm_set1_pd(TWOPOW53_INV_DOUBLE));
    rs = _mm_add_pd(rs, _mm_set1_pd(TWOPOW53_INV_DOUBLE/2.0));
#endif

    return rs;
}

QUALIFIERS void transpose128(__m128i & R0, __m128i & R1, __m128i & R2, __m128i & R3)
{
    __m128i T0, T1, T2, T3;
    T0  = _mm_unpacklo_epi32(R0, R1);
    T1  = _mm_unpacklo_epi32(R2, R3);
    T2  = _mm_unpackhi_epi32(R0, R1);
    T3  = _mm_unpackhi_epi32(R2, R3);
    R0  = _mm_unpacklo_epi64(T0, T1);
    R1  = _mm_unpackhi_epi64(T0, T1);
    R2  = _mm_unpacklo_epi64(T2, T3);
    R3  = _mm_unpackhi_epi64(T2, T3);
}


QUALIFIERS void aesni_float4(__m128i ctr0, __m128i ctr1, __m128i ctr2, __m128i ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m128 & rnd1, __m128 & rnd2, __m128 & rnd3, __m128 & rnd4)
{
    // pack input and call AES
    __m128i k128 = _mm_set_epi32(key3, key2, key1, key0);
    __m128i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);
    for (int i = 0; i < 4; ++i)
    {
        ctr[i] = aesni1xm128i(ctr[i], k128);
    }
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);

    // convert uint32 to float
    rnd1 = _my_cvtepu32_ps(ctr[0]);
    rnd2 = _my_cvtepu32_ps(ctr[1]);
    rnd3 = _my_cvtepu32_ps(ctr[2]);
    rnd4 = _my_cvtepu32_ps(ctr[3]);
    // calculate rnd * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f)
#ifdef __FMA__
    rnd1 = _mm_fmadd_ps(rnd1, _mm_set1_ps(TWOPOW32_INV_FLOAT), _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd2 = _mm_fmadd_ps(rnd2, _mm_set1_ps(TWOPOW32_INV_FLOAT), _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd3 = _mm_fmadd_ps(rnd3, _mm_set1_ps(TWOPOW32_INV_FLOAT), _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd4 = _mm_fmadd_ps(rnd4, _mm_set1_ps(TWOPOW32_INV_FLOAT), _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0));
#else
    rnd1 = _mm_mul_ps(rnd1, _mm_set1_ps(TWOPOW32_INV_FLOAT));
    rnd1 = _mm_add_ps(rnd1, _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
    rnd2 = _mm_mul_ps(rnd2, _mm_set1_ps(TWOPOW32_INV_FLOAT));
    rnd2 = _mm_add_ps(rnd2, _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
    rnd3 = _mm_mul_ps(rnd3, _mm_set1_ps(TWOPOW32_INV_FLOAT));
    rnd3 = _mm_add_ps(rnd3, _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
    rnd4 = _mm_mul_ps(rnd4, _mm_set1_ps(TWOPOW32_INV_FLOAT));
    rnd4 = _mm_add_ps(rnd4, _mm_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
#endif
}


QUALIFIERS void aesni_double2(__m128i ctr0, __m128i ctr1, __m128i ctr2, __m128i ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m128d & rnd1lo, __m128d & rnd1hi, __m128d & rnd2lo, __m128d & rnd2hi)
{
    // pack input and call AES
    __m128i k128 = _mm_set_epi32(key3, key2, key1, key0);
    __m128i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);
    for (int i = 0; i < 4; ++i)
    {
        ctr[i] = aesni1xm128i(ctr[i], k128);
    }
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);

    rnd1lo = _uniform_double_hq<false>(ctr[0], ctr[1]);
    rnd1hi = _uniform_double_hq<true>(ctr[0], ctr[1]);
    rnd2lo = _uniform_double_hq<false>(ctr[2], ctr[3]);
    rnd2hi = _uniform_double_hq<true>(ctr[2], ctr[3]);
}

QUALIFIERS void aesni_float4(uint32 ctr0, __m128i ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             __m128 & rnd1, __m128 & rnd2, __m128 & rnd3, __m128 & rnd4)
{
    __m128i ctr0v = _mm_set1_epi32(ctr0);
    __m128i ctr2v = _mm_set1_epi32(ctr2);
    __m128i ctr3v = _mm_set1_epi32(ctr3);

    aesni_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void aesni_double2(uint32 ctr0, __m128i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m128d & rnd1lo, __m128d & rnd1hi, __m128d & rnd2lo, __m128d & rnd2hi)
{
    __m128i ctr0v = _mm_set1_epi32(ctr0);
    __m128i ctr2v = _mm_set1_epi32(ctr2);
    __m128i ctr3v = _mm_set1_epi32(ctr3);

   aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void aesni_double2(uint32 ctr0, __m128i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m128d & rnd1, __m128d & rnd2)
{
    __m128i ctr0v = _mm_set1_epi32(ctr0);
    __m128i ctr2v = _mm_set1_epi32(ctr2);
    __m128i ctr3v = _mm_set1_epi32(ctr3);

    __m128d ignore;
   aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, ignore, rnd2, ignore);
}
#endif


#ifdef __AVX2__
QUALIFIERS const std::array<__m256i,11> & aesni_roundkeys(const __m256i & k256) {
    alignas(32) std::array<uint32,8> a;
    _mm256_store_si256((__m256i*) a.data(), k256);
    
    static AlignedMap<std::array<uint32,8>, std::array<__m256i,11>> roundkeys;
    
    if(roundkeys.find(a) == roundkeys.end()) {
        auto rk1 = aesni_keygen(_mm256_extractf128_si256(k256, 0));
        auto rk2 = aesni_keygen(_mm256_extractf128_si256(k256, 1));
        for(int i = 0; i < 11; ++i) {
            roundkeys[a][i] = _my256_set_m128i(rk2[i], rk1[i]);
        }
    }
    return roundkeys[a];
}

QUALIFIERS __m256i aesni1xm128i(const __m256i & in, const __m256i & k0) {
#if defined(__VAES__)
    auto k = aesni_roundkeys(k0);
    __m256i x = _mm256_xor_si256(k[0], in);
    x = _mm256_aesenc_epi128(x, k[1]);
    x = _mm256_aesenc_epi128(x, k[2]);
    x = _mm256_aesenc_epi128(x, k[3]);
    x = _mm256_aesenc_epi128(x, k[4]);
    x = _mm256_aesenc_epi128(x, k[5]);
    x = _mm256_aesenc_epi128(x, k[6]);
    x = _mm256_aesenc_epi128(x, k[7]);
    x = _mm256_aesenc_epi128(x, k[8]);
    x = _mm256_aesenc_epi128(x, k[9]);
    x = _mm256_aesenclast_epi128(x, k[10]);
#else
    __m128i a = aesni1xm128i(_mm256_extractf128_si256(in, 0), _mm256_extractf128_si256(k0, 0));
    __m128i b = aesni1xm128i(_mm256_extractf128_si256(in, 1), _mm256_extractf128_si256(k0, 1));
    __m256i x = _my256_set_m128i(b, a);
#endif
    return x;
}

template<bool high>
QUALIFIERS __m256d _uniform_double_hq(__m256i x, __m256i y)
{
    // convert 32 to 64 bit
    if (high)
    {
        x = _mm256_cvtepu32_epi64(_mm256_extracti128_si256(x, 1));
        y = _mm256_cvtepu32_epi64(_mm256_extracti128_si256(y, 1));
    }
    else
    {
        x = _mm256_cvtepu32_epi64(_mm256_extracti128_si256(x, 0));
        y = _mm256_cvtepu32_epi64(_mm256_extracti128_si256(y, 0));
    }

    // calculate z = x ^ y << (53 - 32))
    __m256i z = _mm256_sll_epi64(y, _mm_set1_epi64x(53 - 32));
    z = _mm256_xor_si256(x, z);

    // convert uint64 to double
    __m256d rs = _my256_cvtepu64_pd(z);
    // calculate rs * TWOPOW53_INV_DOUBLE + (TWOPOW53_INV_DOUBLE/2.0)
#ifdef __FMA__
    rs = _mm256_fmadd_pd(rs, _mm256_set1_pd(TWOPOW53_INV_DOUBLE), _mm256_set1_pd(TWOPOW53_INV_DOUBLE/2.0));
#else
    rs = _mm256_mul_pd(rs, _mm256_set1_pd(TWOPOW53_INV_DOUBLE));
    rs = _mm256_add_pd(rs, _mm256_set1_pd(TWOPOW53_INV_DOUBLE/2.0));
#endif

    return rs;
}


QUALIFIERS void transpose128(__m256i & R0, __m256i & R1, __m256i & R2, __m256i & R3)
{
    __m256i T0, T1, T2, T3;
    T0  = _mm256_unpacklo_epi32(R0, R1);
    T1  = _mm256_unpacklo_epi32(R2, R3);
    T2  = _mm256_unpackhi_epi32(R0, R1);
    T3  = _mm256_unpackhi_epi32(R2, R3);
    R0  = _mm256_unpacklo_epi64(T0, T1);
    R1  = _mm256_unpackhi_epi64(T0, T1);
    R2  = _mm256_unpacklo_epi64(T2, T3);
    R3  = _mm256_unpackhi_epi64(T2, T3);
}


QUALIFIERS void aesni_float4(__m256i ctr0, __m256i ctr1, __m256i ctr2, __m256i ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m256 & rnd1, __m256 & rnd2, __m256 & rnd3, __m256 & rnd4)
{
    // pack input and call AES
    __m256i k256 = _mm256_set_epi32(key3, key2, key1, key0, key3, key2, key1, key0);
    __m256i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);
    for (int i = 0; i < 4; ++i)
    {
        ctr[i] = aesni1xm128i(ctr[i], k256);
    }
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);

    // convert uint32 to float
    rnd1 = _my256_cvtepu32_ps(ctr[0]);
    rnd2 = _my256_cvtepu32_ps(ctr[1]);
    rnd3 = _my256_cvtepu32_ps(ctr[2]);
    rnd4 = _my256_cvtepu32_ps(ctr[3]);
    // calculate rnd * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f)
#ifdef __FMA__
    rnd1 = _mm256_fmadd_ps(rnd1, _mm256_set1_ps(TWOPOW32_INV_FLOAT), _mm256_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd2 = _mm256_fmadd_ps(rnd2, _mm256_set1_ps(TWOPOW32_INV_FLOAT), _mm256_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd3 = _mm256_fmadd_ps(rnd3, _mm256_set1_ps(TWOPOW32_INV_FLOAT), _mm256_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd4 = _mm256_fmadd_ps(rnd4, _mm256_set1_ps(TWOPOW32_INV_FLOAT), _mm256_set1_ps(TWOPOW32_INV_FLOAT/2.0));
#else
    rnd1 = _mm256_mul_ps(rnd1, _mm256_set1_ps(TWOPOW32_INV_FLOAT));
    rnd1 = _mm256_add_ps(rnd1, _mm256_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
    rnd2 = _mm256_mul_ps(rnd2, _mm256_set1_ps(TWOPOW32_INV_FLOAT));
    rnd2 = _mm256_add_ps(rnd2, _mm256_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
    rnd3 = _mm256_mul_ps(rnd3, _mm256_set1_ps(TWOPOW32_INV_FLOAT));
    rnd3 = _mm256_add_ps(rnd3, _mm256_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
    rnd4 = _mm256_mul_ps(rnd4, _mm256_set1_ps(TWOPOW32_INV_FLOAT));
    rnd4 = _mm256_add_ps(rnd4, _mm256_set1_ps(TWOPOW32_INV_FLOAT/2.0f));
#endif
}


QUALIFIERS void aesni_double2(__m256i ctr0, __m256i ctr1, __m256i ctr2, __m256i ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m256d & rnd1lo, __m256d & rnd1hi, __m256d & rnd2lo, __m256d & rnd2hi)
{
    // pack input and call AES
    __m256i k256 = _mm256_set_epi32(key3, key2, key1, key0, key3, key2, key1, key0);
    __m256i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);
    for (int i = 0; i < 4; ++i)
    {
        ctr[i] = aesni1xm128i(ctr[i], k256);
    }
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);

    rnd1lo = _uniform_double_hq<false>(ctr[0], ctr[1]);
    rnd1hi = _uniform_double_hq<true>(ctr[0], ctr[1]);
    rnd2lo = _uniform_double_hq<false>(ctr[2], ctr[3]);
    rnd2hi = _uniform_double_hq<true>(ctr[2], ctr[3]);
}

QUALIFIERS void aesni_float4(uint32 ctr0, __m256i ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             __m256 & rnd1, __m256 & rnd2, __m256 & rnd3, __m256 & rnd4)
{
    __m256i ctr0v = _mm256_set1_epi32(ctr0);
    __m256i ctr2v = _mm256_set1_epi32(ctr2);
    __m256i ctr3v = _mm256_set1_epi32(ctr3);

    aesni_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void aesni_double2(uint32 ctr0, __m256i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m256d & rnd1lo, __m256d & rnd1hi, __m256d & rnd2lo, __m256d & rnd2hi)
{
    __m256i ctr0v = _mm256_set1_epi32(ctr0);
    __m256i ctr2v = _mm256_set1_epi32(ctr2);
    __m256i ctr3v = _mm256_set1_epi32(ctr3);

    aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void aesni_double2(uint32 ctr0, __m256i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m256d & rnd1, __m256d & rnd2)
{
#if 0
    __m256i ctr0v = _mm256_set1_epi32(ctr0);
    __m256i ctr2v = _mm256_set1_epi32(ctr2);
    __m256i ctr3v = _mm256_set1_epi32(ctr3);

    __m256d ignore;
    aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, ignore, rnd2, ignore);
#else
    __m128d rnd1lo, rnd1hi, rnd2lo, rnd2hi;
    aesni_double2(ctr0, _mm256_extractf128_si256(ctr1, 0), ctr2, ctr3, key0, key1, key2, key3, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
    rnd1 = _my256_set_m128d(rnd1hi, rnd1lo);
    rnd2 = _my256_set_m128d(rnd2hi, rnd2lo);
#endif
}
#endif


#if defined(__AVX512F__) || defined(__AVX10_512BIT__)
QUALIFIERS const std::array<__m512i,11> & aesni_roundkeys(const __m512i & k512) {
    alignas(64) std::array<uint32,16> a;
    _mm512_store_si512((__m512i*) a.data(), k512);
    
    static AlignedMap<std::array<uint32,16>, std::array<__m512i,11>> roundkeys;
    
    if(roundkeys.find(a) == roundkeys.end()) {
        auto rk1 = aesni_keygen(_mm512_extracti32x4_epi32(k512, 0));
        auto rk2 = aesni_keygen(_mm512_extracti32x4_epi32(k512, 1));
        auto rk3 = aesni_keygen(_mm512_extracti32x4_epi32(k512, 2));
        auto rk4 = aesni_keygen(_mm512_extracti32x4_epi32(k512, 3));
        for(int i = 0; i < 11; ++i) {
            roundkeys[a][i] = _my512_set_m128i(rk4[i], rk3[i], rk2[i], rk1[i]);
        }
    }
    return roundkeys[a];
}

QUALIFIERS __m512i aesni1xm128i(const __m512i & in, const __m512i & k0) {
#ifdef __VAES__
    auto k = aesni_roundkeys(k0);
    __m512i x = _mm512_xor_si512(k[0], in);
    x = _mm512_aesenc_epi128(x, k[1]);
    x = _mm512_aesenc_epi128(x, k[2]);
    x = _mm512_aesenc_epi128(x, k[3]);
    x = _mm512_aesenc_epi128(x, k[4]);
    x = _mm512_aesenc_epi128(x, k[5]);
    x = _mm512_aesenc_epi128(x, k[6]);
    x = _mm512_aesenc_epi128(x, k[7]);
    x = _mm512_aesenc_epi128(x, k[8]);
    x = _mm512_aesenc_epi128(x, k[9]);
    x = _mm512_aesenclast_epi128(x, k[10]);
#else
    __m128i a = aesni1xm128i(_mm512_extracti32x4_epi32(in, 0), _mm512_extracti32x4_epi32(k0, 0));
    __m128i b = aesni1xm128i(_mm512_extracti32x4_epi32(in, 1), _mm512_extracti32x4_epi32(k0, 1));
    __m128i c = aesni1xm128i(_mm512_extracti32x4_epi32(in, 2), _mm512_extracti32x4_epi32(k0, 2));
    __m128i d = aesni1xm128i(_mm512_extracti32x4_epi32(in, 3), _mm512_extracti32x4_epi32(k0, 3));
    __m512i x = _my512_set_m128i(d, c, b, a);
#endif
    return x;
}

template<bool high>
QUALIFIERS __m512d _uniform_double_hq(__m512i x, __m512i y)
{
    // convert 32 to 64 bit
    if (high)
    {
        x = _mm512_cvtepu32_epi64(_mm512_extracti64x4_epi64(x, 1));
        y = _mm512_cvtepu32_epi64(_mm512_extracti64x4_epi64(y, 1));
    }
    else
    {
        x = _mm512_cvtepu32_epi64(_mm512_extracti64x4_epi64(x, 0));
        y = _mm512_cvtepu32_epi64(_mm512_extracti64x4_epi64(y, 0));
    }

    // calculate z = x ^ y << (53 - 32))
    __m512i z = _mm512_sll_epi64(y, _mm_set1_epi64x(53 - 32));
    z = _mm512_xor_si512(x, z);

    // convert uint64 to double
    __m512d rs = _mm512_cvtepu64_pd(z);
    // calculate rs * TWOPOW53_INV_DOUBLE + (TWOPOW53_INV_DOUBLE/2.0)
    rs = _mm512_fmadd_pd(rs, _mm512_set1_pd(TWOPOW53_INV_DOUBLE), _mm512_set1_pd(TWOPOW53_INV_DOUBLE/2.0));

    return rs;
}


QUALIFIERS void transpose128(__m512i & R0, __m512i & R1, __m512i & R2, __m512i & R3)
{
    __m512i T0, T1, T2, T3;
    T0  = _mm512_unpacklo_epi32(R0, R1);
    T1  = _mm512_unpacklo_epi32(R2, R3);
    T2  = _mm512_unpackhi_epi32(R0, R1);
    T3  = _mm512_unpackhi_epi32(R2, R3);
    R0  = _mm512_unpacklo_epi64(T0, T1);
    R1  = _mm512_unpackhi_epi64(T0, T1);
    R2  = _mm512_unpacklo_epi64(T2, T3);
    R3  = _mm512_unpackhi_epi64(T2, T3);
}


QUALIFIERS void aesni_float4(__m512i ctr0, __m512i ctr1, __m512i ctr2, __m512i ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m512 & rnd1, __m512 & rnd2, __m512 & rnd3, __m512 & rnd4)
{
    // pack input and call AES
    __m512i k512 = _mm512_set_epi32(key3, key2, key1, key0, key3, key2, key1, key0,
                                    key3, key2, key1, key0, key3, key2, key1, key0);
    __m512i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);
    for (int i = 0; i < 4; ++i)
    {
        ctr[i] = aesni1xm128i(ctr[i], k512);
    }
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);

    // convert uint32 to float
    rnd1 = _mm512_cvtepu32_ps(ctr[0]);
    rnd2 = _mm512_cvtepu32_ps(ctr[1]);
    rnd3 = _mm512_cvtepu32_ps(ctr[2]);
    rnd4 = _mm512_cvtepu32_ps(ctr[3]);
    // calculate rnd * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f)
    rnd1 = _mm512_fmadd_ps(rnd1, _mm512_set1_ps(TWOPOW32_INV_FLOAT), _mm512_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd2 = _mm512_fmadd_ps(rnd2, _mm512_set1_ps(TWOPOW32_INV_FLOAT), _mm512_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd3 = _mm512_fmadd_ps(rnd3, _mm512_set1_ps(TWOPOW32_INV_FLOAT), _mm512_set1_ps(TWOPOW32_INV_FLOAT/2.0));
    rnd4 = _mm512_fmadd_ps(rnd4, _mm512_set1_ps(TWOPOW32_INV_FLOAT), _mm512_set1_ps(TWOPOW32_INV_FLOAT/2.0));
}


QUALIFIERS void aesni_double2(__m512i ctr0, __m512i ctr1, __m512i ctr2, __m512i ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m512d & rnd1lo, __m512d & rnd1hi, __m512d & rnd2lo, __m512d & rnd2hi)
{
    // pack input and call AES
    __m512i k512 = _mm512_set_epi32(key3, key2, key1, key0, key3, key2, key1, key0,
                                    key3, key2, key1, key0, key3, key2, key1, key0);
    __m512i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);
    for (int i = 0; i < 4; ++i)
    {
        ctr[i] = aesni1xm128i(ctr[i], k512);
    }
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);

    rnd1lo = _uniform_double_hq<false>(ctr[0], ctr[1]);
    rnd1hi = _uniform_double_hq<true>(ctr[0], ctr[1]);
    rnd2lo = _uniform_double_hq<false>(ctr[2], ctr[3]);
    rnd2hi = _uniform_double_hq<true>(ctr[2], ctr[3]);
}

QUALIFIERS void aesni_float4(uint32 ctr0, __m512i ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             __m512 & rnd1, __m512 & rnd2, __m512 & rnd3, __m512 & rnd4)
{
    __m512i ctr0v = _mm512_set1_epi32(ctr0);
    __m512i ctr2v = _mm512_set1_epi32(ctr2);
    __m512i ctr3v = _mm512_set1_epi32(ctr3);

    aesni_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void aesni_double2(uint32 ctr0, __m512i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m512d & rnd1lo, __m512d & rnd1hi, __m512d & rnd2lo, __m512d & rnd2hi)
{
    __m512i ctr0v = _mm512_set1_epi32(ctr0);
    __m512i ctr2v = _mm512_set1_epi32(ctr2);
    __m512i ctr3v = _mm512_set1_epi32(ctr3);

    aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void aesni_double2(uint32 ctr0, __m512i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              __m512d & rnd1, __m512d & rnd2)
{
#if 0
    __m512i ctr0v = _mm512_set1_epi32(ctr0);
    __m512i ctr2v = _mm512_set1_epi32(ctr2);
    __m512i ctr3v = _mm512_set1_epi32(ctr3);

    __m512d ignore;
    aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, ignore, rnd2, ignore);
#else
   __m256d rnd1lo, rnd1hi, rnd2lo, rnd2hi;
   aesni_double2(ctr0, _mm512_extracti64x4_epi64(ctr1, 0), ctr2, ctr3, key0, key1, key2, key3, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
   rnd1 = _my512_set_m256d(rnd1hi, rnd1lo);
   rnd2 = _my512_set_m256d(rnd2hi, rnd2lo);
#endif
}
#endif

#if defined(__ARM_NEON)
QUALIFIERS uint32x4_t aesni_keygen_assist(uint32x4_t temp1, uint32x4_t temp2) {
    uint32x4_t temp3; 
    temp2 = vdupq_laneq_u32(temp2, 3); 
    temp3 = vextq_u32(vdupq_n_u32(0), temp1,  4 - 1);
    temp1 = veorq_u32(temp1, temp3);
    temp3 = vextq_u32(vdupq_n_u32(0), temp3,  4 - 1);
    temp1 = veorq_u32(temp1, temp3);
    temp3 = vextq_u32(vdupq_n_u32(0), temp3,  4 - 1);
    temp1 = veorq_u32(temp1, temp3);
    temp1 = veorq_u32(temp1, temp2); 
    return temp1; 
}

QUALIFIERS uint32x4_t aesni_keygen_assist(uint32x4_t k, const unsigned char imm8)
{
    uint8x16_t a = vreinterpretq_u8_u32(k);
    a = vaeseq_u8(a, vdupq_n_u8(0));
    uint8x16_t dest = {
        a[0x4], a[0x1], a[0xE], a[0xB],
        a[0x1], a[0xE], a[0xB], a[0x4],
        a[0xC], a[0x9], a[0x6], a[0x3],
        a[0x9], a[0x6], a[0x3], a[0xC],
    };
    return vreinterpretq_u32_u8(veorq_u8(dest, vreinterpretq_u8_u32(uint32x4_t{0, imm8, 0, imm8})));
}

QUALIFIERS std::array<uint32x4_t,11> aesni_keygen(uint32x4_t k) {
    std::array<uint32x4_t,11> rk;
    uint32x4_t tmp;
    
    rk[0] = k;
    
    tmp = aesni_keygen_assist(k, 0x1);
    k = aesni_keygen_assist(k, tmp);
    rk[1] = k;
    
    tmp = aesni_keygen_assist(k, 0x2);
    k = aesni_keygen_assist(k, tmp);
    rk[2] = k;
    
    tmp = aesni_keygen_assist(k, 0x4);
    k = aesni_keygen_assist(k, tmp);
    rk[3] = k;
    
    tmp = aesni_keygen_assist(k, 0x8);
    k = aesni_keygen_assist(k, tmp);
    rk[4] = k;
    
    tmp = aesni_keygen_assist(k, 0x10);
    k = aesni_keygen_assist(k, tmp);
    rk[5] = k;
    
    tmp = aesni_keygen_assist(k, 0x20);
    k = aesni_keygen_assist(k, tmp);
    rk[6] = k;
    
    tmp = aesni_keygen_assist(k, 0x40);
    k = aesni_keygen_assist(k, tmp);
    rk[7] = k;
    
    tmp = aesni_keygen_assist(k, 0x80);
    k = aesni_keygen_assist(k, tmp);
    rk[8] = k;
    
    tmp = aesni_keygen_assist(k, 0x1b);
    k = aesni_keygen_assist(k, tmp);
    rk[9] = k;
    
    tmp = aesni_keygen_assist(k, 0x36);
    k = aesni_keygen_assist(k, tmp);
    rk[10] = k;
    
    return rk;
}

QUALIFIERS const std::array<uint32x4_t,11> & aesni_roundkeys(const uint32x4_t & k128) {
    alignas(16) std::array<uint32,4> a;
    vst1q_u32((uint32_t*) a.data(), k128);
    
    static AlignedMap<std::array<uint32,4>, std::array<uint32x4_t,11>> roundkeys;
    
    if(roundkeys.find(a) == roundkeys.end()) {
        auto rk = aesni_keygen(k128);
        roundkeys[a] = rk;
    }
    return roundkeys[a];
}

QUALIFIERS uint32x4_t aesni1xm128i(const uint32x4_t & in, const uint32x4_t & k0) {
    auto k = aesni_roundkeys(k0);
    uint8x16_t x = vaesmcq_u8(vaeseq_u8(vreinterpretq_u8_u32(in), vreinterpretq_u8_u32(k[0])));
    x = vaesmcq_u8(vaeseq_u8(x, vreinterpretq_u8_u32(k[1])));
    x = vaesmcq_u8(vaeseq_u8(x, vreinterpretq_u8_u32(k[2])));
    x = vaesmcq_u8(vaeseq_u8(x, vreinterpretq_u8_u32(k[3])));
    x = vaesmcq_u8(vaeseq_u8(x, vreinterpretq_u8_u32(k[4])));
    x = vaesmcq_u8(vaeseq_u8(x, vreinterpretq_u8_u32(k[5])));
    x = vaesmcq_u8(vaeseq_u8(x, vreinterpretq_u8_u32(k[6])));
    x = vaesmcq_u8(vaeseq_u8(x, vreinterpretq_u8_u32(k[7])));
    x = vaesmcq_u8(vaeseq_u8(x, vreinterpretq_u8_u32(k[8])));
    x = vaeseq_u8(x, vreinterpretq_u8_u32(k[9]));
    x = veorq_u8(x, vreinterpretq_u8_u32(k[10]));
    return vreinterpretq_u32_u8(x);
}

QUALIFIERS void aesni_float4(uint32 ctr0, uint32 ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             float & rnd1, float & rnd2, float & rnd3, float & rnd4)
{
    // pack input and call AES
    uint32x4_t c128 = {ctr0, ctr1, ctr2, ctr3};
    uint32x4_t k128 = {key0, key1, key2, key3};
    c128 = aesni1xm128i(c128, k128);

    // convert uint32 to float
    float32x4_t rs = vcvtq_f32_u32(c128);
    // calculate rs * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f)
    rs = vfmaq_f32(vdupq_n_f32(TWOPOW32_INV_FLOAT/2.0), vdupq_n_f32(TWOPOW32_INV_FLOAT), rs);

    // store result
    rnd1 = vgetq_lane_f32(rs, 0);
    rnd2 = vgetq_lane_f32(rs, 1);
    rnd3 = vgetq_lane_f32(rs, 2);
    rnd4 = vgetq_lane_f32(rs, 3);
}

QUALIFIERS void aesni_double2(uint32 ctr0, uint32 ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              double & rnd1, double & rnd2)
{
    // pack input and call AES
    uint32x4_t c128 = {ctr0, ctr1, ctr2, ctr3};
    uint32x4_t k128 = {key0, key1, key2, key3};
    c128 = aesni1xm128i(c128, k128);

    // convert 32 to 64 bit and put 0th and 2nd element into x, 1st and 3rd element into y
    uint32x4_t x = vandq_u32(c128, uint32x4_t{0xffffffff, 0, 0xffffffff, 0});
    uint32x4_t y = vandq_u32(c128, uint32x4_t{0, 0xffffffff, 0, 0xffffffff});
    y = vextq_u32(y, vdupq_n_u32(0), 1);

    // calculate z = x ^ y << (53 - 32))
    uint64x2_t z = vshlq_n_u64(vreinterpretq_u64_u32(y), 53 - 32);
    z = veorq_u64(vreinterpretq_u64_u32(x), z);

    // convert uint64 to double
    float64x2_t rs = vcvtq_f64_u64(z);
    // calculate rs * TWOPOW53_INV_DOUBLE + (TWOPOW53_INV_DOUBLE/2.0)
    rs = vfmaq_f64(vdupq_n_f64(TWOPOW53_INV_DOUBLE/2.0), vdupq_n_f64(TWOPOW53_INV_DOUBLE), rs);

    // store result
   rnd1 = vgetq_lane_f64(rs, 0);
   rnd2 = vgetq_lane_f64(rs, 1);
}

template<bool high>
QUALIFIERS float64x2_t _uniform_double_hq(uint32x4_t x, uint32x4_t y)
{
    // convert 32 to 64 bit
    if (high)
    {
        x = vzip2q_u32(x, vdupq_n_u32(0));
        y = vzip2q_u32(y, vdupq_n_u32(0));
    }
    else
    {
        x = vzip1q_u32(x, vdupq_n_u32(0));
        y = vzip1q_u32(y, vdupq_n_u32(0));
    }

    // calculate z = x ^ y << (53 - 32))
    uint64x2_t z = vshlq_n_u64(vreinterpretq_u64_u32(y), 53 - 32);
    z = veorq_u64(vreinterpretq_u64_u32(x), z);

    // convert uint64 to double
    float64x2_t rs = vcvtq_f64_u64(z);
    // calculate rs * TWOPOW53_INV_DOUBLE + (TWOPOW53_INV_DOUBLE/2.0)
    rs = vfmaq_f64(vdupq_n_f64(TWOPOW53_INV_DOUBLE/2.0), vdupq_n_f64(TWOPOW53_INV_DOUBLE), rs);

    return rs;
}


QUALIFIERS void transpose128(uint32x4_t & R0, uint32x4_t & R1, uint32x4_t & R2, uint32x4_t & R3)
{
    uint32x4_t T0, T1, T2, T3;
    T0  = vreinterpretq_u32_u64(vtrn1q_u64(vreinterpretq_u64_u32(R0), vreinterpretq_u64_u32(R2)));
    T1  = vreinterpretq_u32_u64(vtrn1q_u64(vreinterpretq_u64_u32(R1), vreinterpretq_u64_u32(R3)));
    T2  = vreinterpretq_u32_u64(vtrn2q_u64(vreinterpretq_u64_u32(R0), vreinterpretq_u64_u32(R2)));
    T3  = vreinterpretq_u32_u64(vtrn2q_u64(vreinterpretq_u64_u32(R1), vreinterpretq_u64_u32(R3)));
    R0  = vtrn1q_u32(T0, T1);
    R1  = vtrn2q_u32(T0, T1);
    R2  = vtrn1q_u32(T2, T3);
    R3  = vtrn2q_u32(T2, T3);
}


QUALIFIERS void aesni_float4(uint32x4_t ctr0, uint32x4_t ctr1, uint32x4_t ctr2, uint32x4_t ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             float32x4_t & rnd1, float32x4_t & rnd2, float32x4_t & rnd3, float32x4_t & rnd4)
{
    uint32x4_t k128 = {key0, key1, key2, key3};
    uint32x4_t ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);
    for (int i = 0; i < 4; ++i)
    {
        ctr[i] = aesni1xm128i(ctr[i], k128);
    }
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);

    // convert uint32 to float
    rnd1 = vcvtq_f32_u32(ctr[0]);
    rnd2 = vcvtq_f32_u32(ctr[1]);
    rnd3 = vcvtq_f32_u32(ctr[2]);
    rnd4 = vcvtq_f32_u32(ctr[3]);
    // calculate rnd * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f)
    rnd1 = vfmaq_f32(vdupq_n_f32(TWOPOW32_INV_FLOAT/2.0), vdupq_n_f32(TWOPOW32_INV_FLOAT), rnd1);
    rnd2 = vfmaq_f32(vdupq_n_f32(TWOPOW32_INV_FLOAT/2.0), vdupq_n_f32(TWOPOW32_INV_FLOAT), rnd2);
    rnd3 = vfmaq_f32(vdupq_n_f32(TWOPOW32_INV_FLOAT/2.0), vdupq_n_f32(TWOPOW32_INV_FLOAT), rnd3);
    rnd4 = vfmaq_f32(vdupq_n_f32(TWOPOW32_INV_FLOAT/2.0), vdupq_n_f32(TWOPOW32_INV_FLOAT), rnd4);
}


QUALIFIERS void aesni_double2(uint32x4_t ctr0, uint32x4_t ctr1, uint32x4_t ctr2, uint32x4_t ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              float64x2_t & rnd1lo, float64x2_t & rnd1hi, float64x2_t & rnd2lo, float64x2_t & rnd2hi)
{
    uint32x4_t k128 = {key0, key1, key2, key3};
    uint32x4_t ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);
    for (int i = 0; i < 4; ++i)
    {
        ctr[i] = aesni1xm128i(ctr[i], k128);
    }
    transpose128(ctr[0], ctr[1], ctr[2], ctr[3]);

    rnd1lo = _uniform_double_hq<false>(ctr[0], ctr[1]);
    rnd1hi = _uniform_double_hq<true>(ctr[0], ctr[1]);
    rnd2lo = _uniform_double_hq<false>(ctr[2], ctr[3]);
    rnd2hi = _uniform_double_hq<true>(ctr[2], ctr[3]);
}

QUALIFIERS void aesni_float4(uint32 ctr0, uint32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             float32x4_t & rnd1, float32x4_t & rnd2, float32x4_t & rnd3, float32x4_t & rnd4)
{
    uint32x4_t ctr0v = vdupq_n_u32(ctr0);
    uint32x4_t ctr2v = vdupq_n_u32(ctr2);
    uint32x4_t ctr3v = vdupq_n_u32(ctr3);

    aesni_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, rnd2, rnd3, rnd4);
}

#ifndef _MSC_VER
QUALIFIERS void aesni_float4(uint32 ctr0, int32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             float32x4_t & rnd1, float32x4_t & rnd2, float32x4_t & rnd3, float32x4_t & rnd4)
{
    aesni_float4(ctr0, vreinterpretq_u32_s32(ctr1), ctr2, ctr3, key0, key1, key2, key3, rnd1, rnd2, rnd3, rnd4);
}
#endif

QUALIFIERS void aesni_double2(uint32 ctr0, uint32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              float64x2_t & rnd1lo, float64x2_t & rnd1hi, float64x2_t & rnd2lo, float64x2_t & rnd2hi)
{
    uint32x4_t ctr0v = vdupq_n_u32(ctr0);
    uint32x4_t ctr2v = vdupq_n_u32(ctr2);
    uint32x4_t ctr3v = vdupq_n_u32(ctr3);

    aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void aesni_double2(uint32 ctr0, uint32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              float64x2_t & rnd1, float64x2_t & rnd2)
{
    uint32x4_t ctr0v = vdupq_n_u32(ctr0);
    uint32x4_t ctr2v = vdupq_n_u32(ctr2);
    uint32x4_t ctr3v = vdupq_n_u32(ctr3);

    float64x2_t ignore;
    aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, ignore, rnd2, ignore);
}

#ifndef _MSC_VER
QUALIFIERS void aesni_double2(uint32 ctr0, int32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              float64x2_t & rnd1, float64x2_t & rnd2)
{
    aesni_double2(ctr0, vreinterpretq_u32_s32(ctr1), ctr2, ctr3, key0, key1, key2, key3, rnd1, rnd2);
}
#endif
#endif

#if defined(__ARM_FEATURE_SVE)
template<bool high>
QUALIFIERS svfloat64_t _uniform_double_hq(svuint32_t x, svuint32_t y)
{
    // convert 32 to 64 bit
    if (high)
    {
        x = svzip2_u32(x, svdup_u32(0));
        y = svzip2_u32(y, svdup_u32(0));
    }
    else
    {
        x = svzip1_u32(x, svdup_u32(0));
        y = svzip1_u32(y, svdup_u32(0));
    }

    // calculate z = x ^ y << (53 - 32))
    svuint64_t z = svlsl_n_u64_x(svptrue_b64(), svreinterpret_u64_u32(y), 53 - 32);
    z = sveor_u64_x(svptrue_b64(), svreinterpret_u64_u32(x), z);

    // convert uint64 to double
    svfloat64_t rs = svcvt_f64_u64_x(svptrue_b64(), z);
    // calculate rs * TWOPOW53_INV_DOUBLE + (TWOPOW53_INV_DOUBLE/2.0)
    rs = svmad_f64_x(svptrue_b64(), rs, svdup_f64(TWOPOW53_INV_DOUBLE), svdup_f64(TWOPOW53_INV_DOUBLE/2.0));

    return rs;
}


QUALIFIERS void transpose128(svuint32x4_t & R)
{
    svuint32_t T0, T1, T2, T3;
    T0  = svreinterpret_u32_u64(svtrn1_u64(svreinterpret_u64_u32(svget4_u32(R, 0)), svreinterpret_u64_u32(svget4_u32(R, 2))));
    T1  = svreinterpret_u32_u64(svtrn1_u64(svreinterpret_u64_u32(svget4_u32(R, 1)), svreinterpret_u64_u32(svget4_u32(R, 3))));
    T2  = svreinterpret_u32_u64(svtrn2_u64(svreinterpret_u64_u32(svget4_u32(R, 0)), svreinterpret_u64_u32(svget4_u32(R, 2))));
    T3  = svreinterpret_u32_u64(svtrn2_u64(svreinterpret_u64_u32(svget4_u32(R, 1)), svreinterpret_u64_u32(svget4_u32(R, 3))));
    R = svset4_u32(R, 0, svtrn1_u32(T0, T1));
    R = svset4_u32(R, 1, svtrn2_u32(T0, T1));
    R = svset4_u32(R, 2, svtrn1_u32(T2, T3));
    R = svset4_u32(R, 3, svtrn2_u32(T2, T3));
}


QUALIFIERS svuint32_t aesni1xm128i(const svuint32_t & in, const uint32x4_t & k0)
{
#ifdef __ARM_FEATURE_SVE2_AES
  auto k = aesni_roundkeys(k0);
  svuint8_t x = svaesmc_u8(svaese_u8(svreinterpret_u8_u32(in), svdup_neonq_u8(vreinterpretq_u8_u32(k[0]))));
  x = svaesmc_u8(svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[1]))));
  x = svaesmc_u8(svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[2]))));
  x = svaesmc_u8(svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[3]))));
  x = svaesmc_u8(svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[4]))));
  x = svaesmc_u8(svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[5]))));
  x = svaesmc_u8(svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[6]))));
  x = svaesmc_u8(svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[7]))));
  x = svaesmc_u8(svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[8]))));
  x = svaese_u8(x, svdup_neonq_u8(vreinterpretq_u8_u32(k[9])));
  x = sveor_u8_x(svptrue_b8(), x, svdup_neonq_u8(vreinterpretq_u8_u32(k[10])));
  return svreinterpret_u32_u8(x);
#else
  svuint32_t x;
  for (int i = 0; i < svcntw(); i += 4)
  {
    svbool_t pred = svbic_z(svptrue_b32(), svwhilelt_b32_u32(0, i+4), svwhilelt_b32_u32(0, i));
    uint32x4_t a = aesni1xm128i(svget_neonq_u32(svcompact_u32(pred, in)), k0);
    x = svsel_u32(pred, svdup_neonq_u32(a), x);
  }
  return x;
#endif
}


QUALIFIERS void aesni_float4(svuint32_t ctr0, svuint32_t ctr1, svuint32_t ctr2, svuint32_t ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             svfloat32_t & rnd1, svfloat32_t & rnd2, svfloat32_t & rnd3, svfloat32_t & rnd4)
{
    uint32x4_t k128 = {key0, key1, key2, key3};
    svuint32x4_t ctr = svcreate4_u32(ctr0, ctr1, ctr2, ctr3);
    transpose128(ctr);
    ctr = svset4_u32(ctr, 0, aesni1xm128i(svget4_u32(ctr, 0), k128));
    ctr = svset4_u32(ctr, 1, aesni1xm128i(svget4_u32(ctr, 1), k128));
    ctr = svset4_u32(ctr, 2, aesni1xm128i(svget4_u32(ctr, 2), k128));
    ctr = svset4_u32(ctr, 3, aesni1xm128i(svget4_u32(ctr, 3), k128));
    transpose128(ctr);

    // convert uint32 to float
    rnd1 = svcvt_f32_u32_x(svptrue_b32(), svget4_u32(ctr, 0));
    rnd2 = svcvt_f32_u32_x(svptrue_b32(), svget4_u32(ctr, 1));
    rnd3 = svcvt_f32_u32_x(svptrue_b32(), svget4_u32(ctr, 2));
    rnd4 = svcvt_f32_u32_x(svptrue_b32(), svget4_u32(ctr, 3));
    // calculate rnd * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f)
    rnd1 = svmad_f32_x(svptrue_b32(), rnd1, svdup_f32(TWOPOW32_INV_FLOAT), svdup_f32(TWOPOW32_INV_FLOAT/2.0));
    rnd2 = svmad_f32_x(svptrue_b32(), rnd2, svdup_f32(TWOPOW32_INV_FLOAT), svdup_f32(TWOPOW32_INV_FLOAT/2.0));
    rnd3 = svmad_f32_x(svptrue_b32(), rnd3, svdup_f32(TWOPOW32_INV_FLOAT), svdup_f32(TWOPOW32_INV_FLOAT/2.0));
    rnd4 = svmad_f32_x(svptrue_b32(), rnd4, svdup_f32(TWOPOW32_INV_FLOAT), svdup_f32(TWOPOW32_INV_FLOAT/2.0));
}


QUALIFIERS void aesni_double2(svuint32_t ctr0, svuint32_t ctr1, svuint32_t ctr2, svuint32_t ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              svfloat64_t & rnd1lo, svfloat64_t & rnd1hi, svfloat64_t & rnd2lo, svfloat64_t & rnd2hi)
{
    uint32x4_t k128 = {key0, key1, key2, key3};
    svuint32x4_t ctr = svcreate4_u32(ctr0, ctr1, ctr2, ctr3);
    transpose128(ctr);
    ctr = svset4_u32(ctr, 0, aesni1xm128i(svget4_u32(ctr, 0), k128));
    ctr = svset4_u32(ctr, 1, aesni1xm128i(svget4_u32(ctr, 1), k128));
    ctr = svset4_u32(ctr, 2, aesni1xm128i(svget4_u32(ctr, 2), k128));
    ctr = svset4_u32(ctr, 3, aesni1xm128i(svget4_u32(ctr, 3), k128));
    transpose128(ctr);

    rnd1lo = _uniform_double_hq<false>(svget4_u32(ctr, 0), svget4_u32(ctr, 1));
    rnd1hi = _uniform_double_hq<true>(svget4_u32(ctr, 0), svget4_u32(ctr, 1));
    rnd2lo = _uniform_double_hq<false>(svget4_u32(ctr, 2), svget4_u32(ctr, 3));
    rnd2hi = _uniform_double_hq<true>(svget4_u32(ctr, 2), svget4_u32(ctr, 3));
}

QUALIFIERS void aesni_float4(uint32 ctr0, svuint32_t ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             svfloat32_t & rnd1, svfloat32_t & rnd2, svfloat32_t & rnd3, svfloat32_t & rnd4)
{
    svuint32_t ctr0v = svdup_u32(ctr0);
    svuint32_t ctr2v = svdup_u32(ctr2);
    svuint32_t ctr3v = svdup_u32(ctr3);

    aesni_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void aesni_float4(uint32 ctr0, svint32_t ctr1, uint32 ctr2, uint32 ctr3,
                             uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                             svfloat32_t & rnd1, svfloat32_t & rnd2, svfloat32_t & rnd3, svfloat32_t & rnd4)
{
    aesni_float4(ctr0, svreinterpret_u32_s32(ctr1), ctr2, ctr3, key0, key1, key2, key3, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void aesni_double2(uint32 ctr0, svuint32_t ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              svfloat64_t & rnd1lo, svfloat64_t & rnd1hi, svfloat64_t & rnd2lo, svfloat64_t & rnd2hi)
{
    svuint32_t ctr0v = svdup_u32(ctr0);
    svuint32_t ctr2v = svdup_u32(ctr2);
    svuint32_t ctr3v = svdup_u32(ctr3);

    aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void aesni_double2(uint32 ctr0, svuint32_t ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              svfloat64_t & rnd1, svfloat64_t & rnd2)
{
    svuint32_t ctr0v = svdup_u32(ctr0);
    svuint32_t ctr2v = svdup_u32(ctr2);
    svuint32_t ctr3v = svdup_u32(ctr3);

    svfloat64_t ignore;
    aesni_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, key2, key3, rnd1, ignore, rnd2, ignore);
}

QUALIFIERS void aesni_double2(uint32 ctr0, svint32_t ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1, uint32 key2, uint32 key3,
                              svfloat64_t & rnd1, svfloat64_t & rnd2)
{
    aesni_double2(ctr0, svreinterpret_u32_s32(ctr1), ctr2, ctr3, key0, key1, key2, key3, rnd1, rnd2);
}
#endif

#undef QUALIFIERS
#undef TWOPOW53_INV_DOUBLE
#undef TWOPOW32_INV_FLOAT
