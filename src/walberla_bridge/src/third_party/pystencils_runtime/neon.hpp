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

#include <array>
#include <cstdint>
#include <arm_neon.h>

#include "./half.h"
#include "./short_arrays.hpp"
#include "./bits/philox_rand.h"

namespace pystencils::runtime::neon
{
    inline int8x8_t vset_s8(int8_t arg0, int8_t arg1, int8_t arg2, int8_t arg3, int8_t arg4, int8_t arg5, int8_t arg6, int8_t arg7)
    {
        alignas(16) std::array<int8_t, 8> values{arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7};
        return vld1_s8(values.data());
    }

    inline int8x16_t vsetq_s8(int8_t arg0, int8_t arg1, int8_t arg2, int8_t arg3, int8_t arg4, int8_t arg5, int8_t arg6, int8_t arg7, int8_t arg8, int8_t arg9, int8_t arg10, int8_t arg11, int8_t arg12, int8_t arg13, int8_t arg14, int8_t arg15)
    {
        alignas(16) std::array<int8_t, 16> values{arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15};
        return vld1q_s8(values.data());
    }

    inline int16x4_t vset_s16(int16_t arg0, int16_t arg1, int16_t arg2, int16_t arg3)
    {
        alignas(16) std::array<int16_t, 4> values{arg0, arg1, arg2, arg3};
        return vld1_s16(values.data());
    }

    inline int16x8_t vsetq_s16(int16_t arg0, int16_t arg1, int16_t arg2, int16_t arg3, int16_t arg4, int16_t arg5, int16_t arg6, int16_t arg7)
    {
        alignas(16) std::array<int16_t, 8> values{arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7};
        return vld1q_s16(values.data());
    }

    inline int32x2_t vset_s32(int32_t arg0, int32_t arg1)
    {
        alignas(16) std::array<int32_t, 2> values{arg0, arg1};
        return vld1_s32(values.data());
    }

    inline int32x4_t vsetq_s32(int32_t arg0, int32_t arg1, int32_t arg2, int32_t arg3)
    {
        alignas(16) std::array<int32_t, 4> values{arg0, arg1, arg2, arg3};
        return vld1q_s32(values.data());
    }

    inline int64x2_t vsetq_s64(int64_t arg0, int64_t arg1)
    {
        alignas(16) std::array<int64_t, 2> values{arg0, arg1};
        return vld1q_s64(values.data());
    }

    inline uint8x8_t vset_u8(uint8_t arg0, uint8_t arg1, uint8_t arg2, uint8_t arg3, uint8_t arg4, uint8_t arg5, uint8_t arg6, uint8_t arg7)
    {
        alignas(16) std::array<uint8_t, 8> values{arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7};
        return vld1_u8(values.data());
    }

    inline uint8x16_t vsetq_u8(uint8_t arg0, uint8_t arg1, uint8_t arg2, uint8_t arg3, uint8_t arg4, uint8_t arg5, uint8_t arg6, uint8_t arg7, uint8_t arg8, uint8_t arg9, uint8_t arg10, uint8_t arg11, uint8_t arg12, uint8_t arg13, uint8_t arg14, uint8_t arg15)
    {
        alignas(16) std::array<uint8_t, 16> values{arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, arg13, arg14, arg15};
        return vld1q_u8(values.data());
    }

    inline uint16x4_t vset_u16(uint16_t arg0, uint16_t arg1, uint16_t arg2, uint16_t arg3)
    {
        alignas(16) std::array<uint16_t, 4> values{arg0, arg1, arg2, arg3};
        return vld1_u16(values.data());
    }

    inline uint16x8_t vsetq_u16(uint16_t arg0, uint16_t arg1, uint16_t arg2, uint16_t arg3, uint16_t arg4, uint16_t arg5, uint16_t arg6, uint16_t arg7)
    {
        alignas(16) std::array<uint16_t, 8> values{arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7};
        return vld1q_u16(values.data());
    }

    inline uint32x2_t vset_u32(uint32_t arg0, uint32_t arg1)
    {
        alignas(16) std::array<uint32_t, 2> values{arg0, arg1};
        return vld1_u32(values.data());
    }

    inline uint32x4_t vsetq_u32(uint32_t arg0, uint32_t arg1, uint32_t arg2, uint32_t arg3)
    {
        alignas(16) std::array<uint32_t, 4> values{arg0, arg1, arg2, arg3};
        return vld1q_u32(values.data());
    }

    inline uint64x2_t vsetq_u64(uint64_t arg0, uint64_t arg1)
    {
        alignas(16) std::array<uint64_t, 2> values{arg0, arg1};
        return vld1q_u64(values.data());
    }

    inline float16x4_t vset_f16(half arg0, half arg1, half arg2, half arg3)
    {
        alignas(16) std::array<half, 4> values{arg0, arg1, arg2, arg3};
        return vld1_f16((const float16_t *)values.data());
    }

    inline float16x8_t vsetq_f16(half arg0, half arg1, half arg2, half arg3, half arg4, half arg5, half arg6, half arg7)
    {
        alignas(16) std::array<half, 8> values{arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7};
        return vld1q_f16((const float16_t *)values.data());
    }

    inline float32x2_t vset_f32(float arg0, float arg1)
    {
        alignas(16) std::array<float, 2> values{arg0, arg1};
        return vld1_f32(values.data());
    }

    inline float32x4_t vsetq_f32(float arg0, float arg1, float arg2, float arg3)
    {
        alignas(16) std::array<float, 4> values{arg0, arg1, arg2, arg3};
        return vld1q_f32(values.data());
    }

    inline float64x2_t vsetq_f64(double arg0, double arg1)
    {
        alignas(16) std::array<double, 2> values{arg0, arg1};
        return vld1q_f64(values.data());
    }

    //  Horizontal multiplication - unsigned integers

    inline uint16_t vmulv_u16(uint16x4_t v)
    {
        return (vget_lane_u16(v, 0) * vget_lane_u16(v, 1)) * (vget_lane_u16(v, 2) * vget_lane_u16(v, 3));
    }

    inline uint16_t vmulvq_u16(uint16x8_t v)
    {
        return ((vgetq_lane_u16(v, 0) * vgetq_lane_u16(v, 1)) * (vgetq_lane_u16(v, 2) * vgetq_lane_u16(v, 3))) * ((vgetq_lane_u16(v, 4) * vgetq_lane_u16(v, 5)) * (vgetq_lane_u16(v, 6) * vgetq_lane_u16(v, 7)));
    }

    inline uint32_t vmulv_u32(uint32x2_t v)
    {
        return vget_lane_u32(v, 0) * vget_lane_u32(v, 1);
    }

    inline uint32_t vmulvq_u32(uint32x4_t v)
    {
        return (vgetq_lane_u32(v, 0) * vgetq_lane_u32(v, 1)) * (vgetq_lane_u32(v, 2) * vgetq_lane_u32(v, 3));
    }

    inline uint64_t vmulvq_u64(uint64x2_t v)
    {
        return vgetq_lane_u64(v, 0) * vgetq_lane_u64(v, 1);
    }

    //  Horizontal multiplication - signed integers

    inline int16_t vmulv_s16(int16x4_t v)
    {
        return (vget_lane_s16(v, 0) * vget_lane_s16(v, 1)) * (vget_lane_s16(v, 2) * vget_lane_s16(v, 3));
    }

    inline int16_t vmulvq_s16(int16x8_t v)
    {
        return ((vgetq_lane_s16(v, 0) * vgetq_lane_s16(v, 1)) * (vgetq_lane_s16(v, 2) * vgetq_lane_s16(v, 3))) * ((vgetq_lane_s16(v, 4) * vgetq_lane_s16(v, 5)) * (vgetq_lane_s16(v, 6) * vgetq_lane_s16(v, 7)));
    }

    inline int32_t vmulv_s32(int32x2_t v)
    {
        return vget_lane_s32(v, 0) * vget_lane_s32(v, 1);
    }

    inline int32_t vmulvq_s32(int32x4_t v)
    {
        return (vgetq_lane_s32(v, 0) * vgetq_lane_s32(v, 1)) * (vgetq_lane_s32(v, 2) * vgetq_lane_s32(v, 3));
    }

    inline int64_t vmulvq_s64(int64x2_t v)
    {
        return vgetq_lane_s64(v, 0) * vgetq_lane_s64(v, 1);
    }

    //  Horizontal multiplication - floating point

    inline half vmulv_f16(float16x4_t v)
    {
        return ((half)(vget_lane_f16(v, 0)) * (half)(vget_lane_f16(v, 1))) * ((half)(vget_lane_f16(v, 2)) * (half)(vget_lane_f16(v, 3)));
    }

    inline half vmulvq_f16(float16x8_t v)
    {
        return (((half)(vgetq_lane_f16(v, 0)) * (half)(vgetq_lane_f16(v, 1))) * ((half)(vgetq_lane_f16(v, 2)) * (half)(vgetq_lane_f16(v, 3)))) * (((half)(vgetq_lane_f16(v, 4)) * (half)(vgetq_lane_f16(v, 5))) * ((half)(vgetq_lane_f16(v, 6)) * (half)(vgetq_lane_f16(v, 7))));
    }

    inline float vmulv_f32(float32x2_t v)
    {
        return vget_lane_f32(v, 0) * vget_lane_f32(v, 1);
    }

    inline float vmulvq_f32(float32x4_t v)
    {
        return (vgetq_lane_f32(v, 0) * vgetq_lane_f32(v, 1)) * (vgetq_lane_f32(v, 2) * vgetq_lane_f32(v, 3));
    }

    inline double vmulvq_f64(float64x2_t v)
    {
        return vgetq_lane_f64(v, 0) * vgetq_lane_f64(v, 1);
    }

    namespace random
    {
        template <typename Idx>
        ShortArray<float32x4_t, 4> philox_float32(uint32x4_t ctr0, uint32x4_t ctr1, uint32x4_t ctr2, uint32x4_t ctr3, Idx key0, Idx key1)
        {
            ShortArray<float32x4_t, 4> result;
            pystencils::runtime::random::detail::philox_float4(
                ctr0, ctr1, ctr2, ctr3,                    //
                (uint32_t)key0, (uint32_t)key1,            //
                result[0], result[1], result[2], result[3] //
            );
            return result;
        }

        template <typename Idx>
        ShortArray<float32x4_t, 4> philox_float32(int32x4_t ctr0, int32x4_t ctr1, int32x4_t ctr2, int32x4_t ctr3, Idx key0, Idx key1)
        {
            return philox_float32(
                vreinterpretq_u32_s32(ctr0),
                vreinterpretq_u32_s32(ctr1),
                vreinterpretq_u32_s32(ctr2),
                vreinterpretq_u32_s32(ctr3),
                key0,
                key1);
        }

        template <typename Idx>
        ShortArray<float64x2_t, 2> philox_float64(uint64x2_t ctr0, uint64x2_t ctr1, uint64x2_t ctr2, uint64x2_t ctr3, Idx key0, Idx key1)
        {
            uint32x4_t ctr0_u32{vuzp1q_u32(vreinterpretq_u32_s64(ctr0), vdupq_n_u32(0u))};
            uint32x4_t ctr1_u32{vuzp1q_u32(vreinterpretq_u32_s64(ctr1), vdupq_n_u32(0u))};
            uint32x4_t ctr2_u32{vuzp1q_u32(vreinterpretq_u32_s64(ctr2), vdupq_n_u32(0u))};
            uint32x4_t ctr3_u32{vuzp1q_u32(vreinterpretq_u32_s64(ctr3), vdupq_n_u32(0u))};

            ShortArray<float64x2_t, 2> result;
            float64x2_t ignore;
            pystencils::runtime::random::detail::philox_double2(
                ctr0_u32, ctr1_u32, ctr2_u32, ctr3_u32, //
                (uint32_t)key0, (uint32_t)key1,         //
                result[0], ignore, result[1], ignore    //
            );
            return result;
        }

        template <typename Idx>
        ShortArray<float64x2_t, 2> philox_float64(int64x2_t ctr0, int64x2_t ctr1, int64x2_t ctr2, int64x2_t ctr3, Idx key0, Idx key1)
        {
            return philox_float64(
                vreinterpretq_u64_s64(ctr0),
                vreinterpretq_u64_s64(ctr1),
                vreinterpretq_u64_s64(ctr2),
                vreinterpretq_u64_s64(ctr3),
                key0,
                key1);
        }
    } // namespace random

} // namespace pystencils::runtime::neon
