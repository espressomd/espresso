#include <cstdint>

#ifndef __CUDA_ARCH__
#define QUALIFIERS inline
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