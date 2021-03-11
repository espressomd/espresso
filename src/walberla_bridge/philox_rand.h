#include <cstdint>

#if defined(__SSE2__) || defined(_MSC_VER)
#include <emmintrin.h> // SSE2
#endif
#ifdef __AVX2__
#include <immintrin.h> // AVX*
#elif defined(__SSE4_1__) || defined(_MSC_VER)
#include <smmintrin.h>  // SSE4
#ifdef __FMA__
#include <immintrin.h> // FMA
#endif
#endif

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

#ifndef __CUDA_ARCH__
#define QUALIFIERS inline
#include "myintrin.h"
#else
#define QUALIFIERS static __forceinline__ __device__
#endif

#define PHILOX_W32_0   (0x9E3779B9)
#define PHILOX_W32_1   (0xBB67AE85)
#define PHILOX_M4x32_0 (0xD2511F53)
#define PHILOX_M4x32_1 (0xCD9E8D57)
#define TWOPOW53_INV_DOUBLE (1.1102230246251565e-16)
#define TWOPOW32_INV_FLOAT (2.3283064e-10f)

typedef std::uint32_t uint32;
typedef std::uint64_t uint64;


QUALIFIERS uint32 mulhilo32(uint32 a, uint32 b, uint32* hip)
{
#ifndef __CUDA_ARCH__
    // host code
    uint64 product = ((uint64)a) * ((uint64)b);
    *hip = product >> 32;
    return (uint32)product;
#else
    // device code
    *hip = __umulhi(a,b);
    return a*b;
#endif
}

QUALIFIERS void _philox4x32round(uint32* ctr, uint32* key)
{
    uint32 hi0;
    uint32 hi1;
    uint32 lo0 = mulhilo32(PHILOX_M4x32_0, ctr[0], &hi0);
    uint32 lo1 = mulhilo32(PHILOX_M4x32_1, ctr[2], &hi1);

    ctr[0] = hi1^ctr[1]^key[0];
    ctr[1] = lo1;
    ctr[2] = hi0^ctr[3]^key[1];
    ctr[3] = lo0;
}

QUALIFIERS void _philox4x32bumpkey(uint32* key)
{
    key[0] += PHILOX_W32_0;
    key[1] += PHILOX_W32_1;
}

QUALIFIERS double _uniform_double_hq(uint32 x, uint32 y)
{
    uint64 z = (uint64)x ^ ((uint64)y << (53 - 32));
    return z * TWOPOW53_INV_DOUBLE + (TWOPOW53_INV_DOUBLE/2.0);
}


QUALIFIERS void philox_double2(uint32 ctr0, uint32 ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1, double & rnd1, double & rnd2)
{
    uint32 key[2] = {key0, key1};
    uint32 ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

    rnd1 = _uniform_double_hq(ctr[0], ctr[1]);
    rnd2 = _uniform_double_hq(ctr[2], ctr[3]);
}


QUALIFIERS void philox_float4(uint32 ctr0, uint32 ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1,
                              float & rnd1, float & rnd2, float & rnd3, float & rnd4)
{
    uint32 key[2] = {key0, key1};
    uint32 ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

    rnd1 = ctr[0] * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f);
    rnd2 = ctr[1] * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f);
    rnd3 = ctr[2] * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f);
    rnd4 = ctr[3] * TWOPOW32_INV_FLOAT + (TWOPOW32_INV_FLOAT/2.0f);
}

#ifndef __CUDA_ARCH__
#if defined(__SSE4_1__) || defined(_MSC_VER)
QUALIFIERS void _philox4x32round(__m128i* ctr, __m128i* key)
{
    __m128i lohi0a = _mm_mul_epu32(ctr[0], _mm_set1_epi32(PHILOX_M4x32_0));
    __m128i lohi0b = _mm_mul_epu32(_mm_srli_epi64(ctr[0], 32), _mm_set1_epi32(PHILOX_M4x32_0));
    __m128i lohi1a = _mm_mul_epu32(ctr[2], _mm_set1_epi32(PHILOX_M4x32_1));
    __m128i lohi1b = _mm_mul_epu32(_mm_srli_epi64(ctr[2], 32), _mm_set1_epi32(PHILOX_M4x32_1));

    lohi0a = _mm_shuffle_epi32(lohi0a, 0xD8);
    lohi0b = _mm_shuffle_epi32(lohi0b, 0xD8);
    lohi1a = _mm_shuffle_epi32(lohi1a, 0xD8);
    lohi1b = _mm_shuffle_epi32(lohi1b, 0xD8);

    __m128i lo0 = _mm_unpacklo_epi32(lohi0a, lohi0b);
    __m128i hi0 = _mm_unpackhi_epi32(lohi0a, lohi0b);
    __m128i lo1 = _mm_unpacklo_epi32(lohi1a, lohi1b);
    __m128i hi1 = _mm_unpackhi_epi32(lohi1a, lohi1b);

    ctr[0] = _mm_xor_si128(_mm_xor_si128(hi1, ctr[1]), key[0]);
    ctr[1] = lo1;
    ctr[2] = _mm_xor_si128(_mm_xor_si128(hi0, ctr[3]), key[1]);
    ctr[3] = lo0;
}

QUALIFIERS void _philox4x32bumpkey(__m128i* key)
{
    key[0] = _mm_add_epi32(key[0], _mm_set1_epi32(PHILOX_W32_0));
    key[1] = _mm_add_epi32(key[1], _mm_set1_epi32(PHILOX_W32_1));
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


QUALIFIERS void philox_float4(__m128i ctr0, __m128i ctr1, __m128i ctr2, __m128i ctr3,
                              uint32 key0, uint32 key1,
                              __m128 & rnd1, __m128 & rnd2, __m128 & rnd3, __m128 & rnd4)
{
    __m128i key[2] = {_mm_set1_epi32(key0), _mm_set1_epi32(key1)};
    __m128i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

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


QUALIFIERS void philox_double2(__m128i ctr0, __m128i ctr1, __m128i ctr2, __m128i ctr3,
                               uint32 key0, uint32 key1,
                               __m128d & rnd1lo, __m128d & rnd1hi, __m128d & rnd2lo, __m128d & rnd2hi)
{
    __m128i key[2] = {_mm_set1_epi32(key0), _mm_set1_epi32(key1)};
    __m128i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

    rnd1lo = _uniform_double_hq<false>(ctr[0], ctr[1]);
    rnd1hi = _uniform_double_hq<true>(ctr[0], ctr[1]);
    rnd2lo = _uniform_double_hq<false>(ctr[2], ctr[3]);
    rnd2hi = _uniform_double_hq<true>(ctr[2], ctr[3]);
}

QUALIFIERS void philox_float4(uint32 ctr0, __m128i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1,
                              __m128 & rnd1, __m128 & rnd2, __m128 & rnd3, __m128 & rnd4)
{
    __m128i ctr0v = _mm_set1_epi32(ctr0);
    __m128i ctr2v = _mm_set1_epi32(ctr2);
    __m128i ctr3v = _mm_set1_epi32(ctr3);

    philox_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void philox_double2(uint32 ctr0, __m128i ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               __m128d & rnd1lo, __m128d & rnd1hi, __m128d & rnd2lo, __m128d & rnd2hi)
{
    __m128i ctr0v = _mm_set1_epi32(ctr0);
    __m128i ctr2v = _mm_set1_epi32(ctr2);
    __m128i ctr3v = _mm_set1_epi32(ctr3);

    philox_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void philox_double2(uint32 ctr0, __m128i ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               __m128d & rnd1, __m128d & rnd2)
{
    __m128i ctr0v = _mm_set1_epi32(ctr0);
    __m128i ctr2v = _mm_set1_epi32(ctr2);
    __m128i ctr3v = _mm_set1_epi32(ctr3);

    __m128d ignore;
    philox_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1, ignore, rnd2, ignore);
}
#endif

#if defined(__ARM_NEON)
QUALIFIERS void _philox4x32round(uint32x4_t* ctr, uint32x4_t* key)
{
    uint32x4_t lohi0a = vreinterpretq_u32_u64(vmull_u32(vget_low_u32(ctr[0]), vdup_n_u32(PHILOX_M4x32_0)));
    uint32x4_t lohi0b = vreinterpretq_u32_u64(vmull_high_u32(ctr[0], vdupq_n_u32(PHILOX_M4x32_0)));
    uint32x4_t lohi1a = vreinterpretq_u32_u64(vmull_u32(vget_low_u32(ctr[2]), vdup_n_u32(PHILOX_M4x32_1)));
    uint32x4_t lohi1b = vreinterpretq_u32_u64(vmull_high_u32(ctr[2], vdupq_n_u32(PHILOX_M4x32_1)));

    uint32x4_t lo0 = vuzp1q_u32(lohi0a, lohi0b);
    uint32x4_t lo1 = vuzp1q_u32(lohi1a, lohi1b);
    uint32x4_t hi0 = vuzp2q_u32(lohi0a, lohi0b);
    uint32x4_t hi1 = vuzp2q_u32(lohi1a, lohi1b);

    ctr[0] = veorq_u32(veorq_u32(hi1, ctr[1]), key[0]);
    ctr[1] = lo1;
    ctr[2] = veorq_u32(veorq_u32(hi0, ctr[3]), key[1]);
    ctr[3] = lo0;
}

QUALIFIERS void _philox4x32bumpkey(uint32x4_t* key)
{
    key[0] = vaddq_u32(key[0], vdupq_n_u32(PHILOX_W32_0));
    key[1] = vaddq_u32(key[1], vdupq_n_u32(PHILOX_W32_1));
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


QUALIFIERS void philox_float4(uint32x4_t ctr0, uint32x4_t ctr1, uint32x4_t ctr2, uint32x4_t ctr3,
                              uint32 key0, uint32 key1,
                              float32x4_t & rnd1, float32x4_t & rnd2, float32x4_t & rnd3, float32x4_t & rnd4)
{
    uint32x4_t key[2] = {vdupq_n_u32(key0), vdupq_n_u32(key1)};
    uint32x4_t ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

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


QUALIFIERS void philox_double2(uint32x4_t ctr0, uint32x4_t ctr1, uint32x4_t ctr2, uint32x4_t ctr3,
                               uint32 key0, uint32 key1,
                               float64x2_t & rnd1lo, float64x2_t & rnd1hi, float64x2_t & rnd2lo, float64x2_t & rnd2hi)
{
    uint32x4_t key[2] = {vdupq_n_u32(key0), vdupq_n_u32(key1)};
    uint32x4_t ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

    rnd1lo = _uniform_double_hq<false>(ctr[0], ctr[1]);
    rnd1hi = _uniform_double_hq<true>(ctr[0], ctr[1]);
    rnd2lo = _uniform_double_hq<false>(ctr[2], ctr[3]);
    rnd2hi = _uniform_double_hq<true>(ctr[2], ctr[3]);
}

QUALIFIERS void philox_float4(uint32 ctr0, uint32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1,
                              float32x4_t & rnd1, float32x4_t & rnd2, float32x4_t & rnd3, float32x4_t & rnd4)
{
    uint32x4_t ctr0v = vdupq_n_u32(ctr0);
    uint32x4_t ctr2v = vdupq_n_u32(ctr2);
    uint32x4_t ctr3v = vdupq_n_u32(ctr3);

    philox_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void philox_float4(uint32 ctr0, int32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1,
                              float32x4_t & rnd1, float32x4_t & rnd2, float32x4_t & rnd3, float32x4_t & rnd4)
{
    philox_float4(ctr0, vreinterpretq_u32_s32(ctr1), ctr2, ctr3, key0, key1, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void philox_double2(uint32 ctr0, uint32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               float64x2_t & rnd1lo, float64x2_t & rnd1hi, float64x2_t & rnd2lo, float64x2_t & rnd2hi)
{
    uint32x4_t ctr0v = vdupq_n_u32(ctr0);
    uint32x4_t ctr2v = vdupq_n_u32(ctr2);
    uint32x4_t ctr3v = vdupq_n_u32(ctr3);

    philox_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void philox_double2(uint32 ctr0, uint32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               float64x2_t & rnd1, float64x2_t & rnd2)
{
    uint32x4_t ctr0v = vdupq_n_u32(ctr0);
    uint32x4_t ctr2v = vdupq_n_u32(ctr2);
    uint32x4_t ctr3v = vdupq_n_u32(ctr3);

    float64x2_t ignore;
    philox_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1, ignore, rnd2, ignore);
}

QUALIFIERS void philox_double2(uint32 ctr0, int32x4_t ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               float64x2_t & rnd1, float64x2_t & rnd2)
{
    philox_double2(ctr0, vreinterpretq_u32_s32(ctr1), ctr2, ctr3, key0, key1, rnd1, rnd2);
}
#endif

#ifdef __AVX2__
QUALIFIERS void _philox4x32round(__m256i* ctr, __m256i* key)
{
    __m256i lohi0a = _mm256_mul_epu32(ctr[0], _mm256_set1_epi32(PHILOX_M4x32_0));
    __m256i lohi0b = _mm256_mul_epu32(_mm256_srli_epi64(ctr[0], 32), _mm256_set1_epi32(PHILOX_M4x32_0));
    __m256i lohi1a = _mm256_mul_epu32(ctr[2], _mm256_set1_epi32(PHILOX_M4x32_1));
    __m256i lohi1b = _mm256_mul_epu32(_mm256_srli_epi64(ctr[2], 32), _mm256_set1_epi32(PHILOX_M4x32_1));

    lohi0a = _mm256_shuffle_epi32(lohi0a, 0xD8);
    lohi0b = _mm256_shuffle_epi32(lohi0b, 0xD8);
    lohi1a = _mm256_shuffle_epi32(lohi1a, 0xD8);
    lohi1b = _mm256_shuffle_epi32(lohi1b, 0xD8);

    __m256i lo0 = _mm256_unpacklo_epi32(lohi0a, lohi0b);
    __m256i hi0 = _mm256_unpackhi_epi32(lohi0a, lohi0b);
    __m256i lo1 = _mm256_unpacklo_epi32(lohi1a, lohi1b);
    __m256i hi1 = _mm256_unpackhi_epi32(lohi1a, lohi1b);

    ctr[0] = _mm256_xor_si256(_mm256_xor_si256(hi1, ctr[1]), key[0]);
    ctr[1] = lo1;
    ctr[2] = _mm256_xor_si256(_mm256_xor_si256(hi0, ctr[3]), key[1]);
    ctr[3] = lo0;
}

QUALIFIERS void _philox4x32bumpkey(__m256i* key)
{
    key[0] = _mm256_add_epi32(key[0], _mm256_set1_epi32(PHILOX_W32_0));
    key[1] = _mm256_add_epi32(key[1], _mm256_set1_epi32(PHILOX_W32_1));
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


QUALIFIERS void philox_float4(__m256i ctr0, __m256i ctr1, __m256i ctr2, __m256i ctr3,
                              uint32 key0, uint32 key1,
                              __m256 & rnd1, __m256 & rnd2, __m256 & rnd3, __m256 & rnd4)
{
    __m256i key[2] = {_mm256_set1_epi32(key0), _mm256_set1_epi32(key1)};
    __m256i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

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


QUALIFIERS void philox_double2(__m256i ctr0, __m256i ctr1, __m256i ctr2, __m256i ctr3,
                               uint32 key0, uint32 key1,
                               __m256d & rnd1lo, __m256d & rnd1hi, __m256d & rnd2lo, __m256d & rnd2hi)
{
    __m256i key[2] = {_mm256_set1_epi32(key0), _mm256_set1_epi32(key1)};
    __m256i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

    rnd1lo = _uniform_double_hq<false>(ctr[0], ctr[1]);
    rnd1hi = _uniform_double_hq<true>(ctr[0], ctr[1]);
    rnd2lo = _uniform_double_hq<false>(ctr[2], ctr[3]);
    rnd2hi = _uniform_double_hq<true>(ctr[2], ctr[3]);
}

QUALIFIERS void philox_float4(uint32 ctr0, __m256i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1,
                              __m256 & rnd1, __m256 & rnd2, __m256 & rnd3, __m256 & rnd4)
{
    __m256i ctr0v = _mm256_set1_epi32(ctr0);
    __m256i ctr2v = _mm256_set1_epi32(ctr2);
    __m256i ctr3v = _mm256_set1_epi32(ctr3);

    philox_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void philox_double2(uint32 ctr0, __m256i ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               __m256d & rnd1lo, __m256d & rnd1hi, __m256d & rnd2lo, __m256d & rnd2hi)
{
    __m256i ctr0v = _mm256_set1_epi32(ctr0);
    __m256i ctr2v = _mm256_set1_epi32(ctr2);
    __m256i ctr3v = _mm256_set1_epi32(ctr3);

    philox_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void philox_double2(uint32 ctr0, __m256i ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               __m256d & rnd1, __m256d & rnd2)
{
#if 0
    __m256i ctr0v = _mm256_set1_epi32(ctr0);
    __m256i ctr2v = _mm256_set1_epi32(ctr2);
    __m256i ctr3v = _mm256_set1_epi32(ctr3);

    __m256d ignore;
    philox_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1, ignore, rnd2, ignore);
#else
    __m128d rnd1lo, rnd1hi, rnd2lo, rnd2hi;
    philox_double2(ctr0, _mm256_extractf128_si256(ctr1, 0), ctr2, ctr3, key0, key1, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
    rnd1 = _my256_set_m128d(rnd1hi, rnd1lo);
    rnd2 = _my256_set_m128d(rnd2hi, rnd2lo);
#endif
}
#endif

#ifdef __AVX512F__
QUALIFIERS void _philox4x32round(__m512i* ctr, __m512i* key)
{
    __m512i lohi0a = _mm512_mul_epu32(ctr[0], _mm512_set1_epi32(PHILOX_M4x32_0));
    __m512i lohi0b = _mm512_mul_epu32(_mm512_srli_epi64(ctr[0], 32), _mm512_set1_epi32(PHILOX_M4x32_0));
    __m512i lohi1a = _mm512_mul_epu32(ctr[2], _mm512_set1_epi32(PHILOX_M4x32_1));
    __m512i lohi1b = _mm512_mul_epu32(_mm512_srli_epi64(ctr[2], 32), _mm512_set1_epi32(PHILOX_M4x32_1));

    lohi0a = _mm512_shuffle_epi32(lohi0a, _MM_PERM_DBCA);
    lohi0b = _mm512_shuffle_epi32(lohi0b, _MM_PERM_DBCA);
    lohi1a = _mm512_shuffle_epi32(lohi1a, _MM_PERM_DBCA);
    lohi1b = _mm512_shuffle_epi32(lohi1b, _MM_PERM_DBCA);

    __m512i lo0 = _mm512_unpacklo_epi32(lohi0a, lohi0b);
    __m512i hi0 = _mm512_unpackhi_epi32(lohi0a, lohi0b);
    __m512i lo1 = _mm512_unpacklo_epi32(lohi1a, lohi1b);
    __m512i hi1 = _mm512_unpackhi_epi32(lohi1a, lohi1b);

    ctr[0] = _mm512_xor_si512(_mm512_xor_si512(hi1, ctr[1]), key[0]);
    ctr[1] = lo1;
    ctr[2] = _mm512_xor_si512(_mm512_xor_si512(hi0, ctr[3]), key[1]);
    ctr[3] = lo0;
}

QUALIFIERS void _philox4x32bumpkey(__m512i* key)
{
    key[0] = _mm512_add_epi32(key[0], _mm512_set1_epi32(PHILOX_W32_0));
    key[1] = _mm512_add_epi32(key[1], _mm512_set1_epi32(PHILOX_W32_1));
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


QUALIFIERS void philox_float4(__m512i ctr0, __m512i ctr1, __m512i ctr2, __m512i ctr3,
                              uint32 key0, uint32 key1,
                              __m512 & rnd1, __m512 & rnd2, __m512 & rnd3, __m512 & rnd4)
{
    __m512i key[2] = {_mm512_set1_epi32(key0), _mm512_set1_epi32(key1)};
    __m512i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

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


QUALIFIERS void philox_double2(__m512i ctr0, __m512i ctr1, __m512i ctr2, __m512i ctr3,
                               uint32 key0, uint32 key1,
                               __m512d & rnd1lo, __m512d & rnd1hi, __m512d & rnd2lo, __m512d & rnd2hi)
{
    __m512i key[2] = {_mm512_set1_epi32(key0), _mm512_set1_epi32(key1)};
    __m512i ctr[4] = {ctr0, ctr1, ctr2, ctr3};
    _philox4x32round(ctr, key);                           // 1
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 2
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 3
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 4
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 5
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 6
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 7
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 8
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 9
    _philox4x32bumpkey(key); _philox4x32round(ctr, key);  // 10

    rnd1lo = _uniform_double_hq<false>(ctr[0], ctr[1]);
    rnd1hi = _uniform_double_hq<true>(ctr[0], ctr[1]);
    rnd2lo = _uniform_double_hq<false>(ctr[2], ctr[3]);
    rnd2hi = _uniform_double_hq<true>(ctr[2], ctr[3]);
}

QUALIFIERS void philox_float4(uint32 ctr0, __m512i ctr1, uint32 ctr2, uint32 ctr3,
                              uint32 key0, uint32 key1,
                              __m512 & rnd1, __m512 & rnd2, __m512 & rnd3, __m512 & rnd4)
{
    __m512i ctr0v = _mm512_set1_epi32(ctr0);
    __m512i ctr2v = _mm512_set1_epi32(ctr2);
    __m512i ctr3v = _mm512_set1_epi32(ctr3);

    philox_float4(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1, rnd2, rnd3, rnd4);
}

QUALIFIERS void philox_double2(uint32 ctr0, __m512i ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               __m512d & rnd1lo, __m512d & rnd1hi, __m512d & rnd2lo, __m512d & rnd2hi)
{
    __m512i ctr0v = _mm512_set1_epi32(ctr0);
    __m512i ctr2v = _mm512_set1_epi32(ctr2);
    __m512i ctr3v = _mm512_set1_epi32(ctr3);

    philox_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
}

QUALIFIERS void philox_double2(uint32 ctr0, __m512i ctr1, uint32 ctr2, uint32 ctr3,
                               uint32 key0, uint32 key1,
                               __m512d & rnd1, __m512d & rnd2)
{
#if 0
    __m512i ctr0v = _mm512_set1_epi32(ctr0);
    __m512i ctr2v = _mm512_set1_epi32(ctr2);
    __m512i ctr3v = _mm512_set1_epi32(ctr3);

    __m512d ignore;
    philox_double2(ctr0v, ctr1, ctr2v, ctr3v, key0, key1, rnd1, ignore, rnd2, ignore);
#else
   __m256d rnd1lo, rnd1hi, rnd2lo, rnd2hi;
   philox_double2(ctr0, _mm512_extracti64x4_epi64(ctr1, 0), ctr2, ctr3, key0, key1, rnd1lo, rnd1hi, rnd2lo, rnd2hi);
   rnd1 = _my512_set_m256d(rnd1hi, rnd1lo);
   rnd2 = _my512_set_m256d(rnd2hi, rnd2lo);
#endif
}
#endif
#endif

