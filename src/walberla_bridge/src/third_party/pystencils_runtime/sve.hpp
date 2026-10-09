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

#include <arm_sve.h>

#include "./half.h"
#include "./bits/philox_rand.h"

namespace pystencils::runtime::sve
{
    namespace random
    {
        template <typename Idx>
        svfloat32x4_t philox_float32(svuint32_t ctr0, svuint32_t ctr1, svuint32_t ctr2, svuint32_t ctr3, Idx key0, Idx key1)
        {
            svfloat32_t result0, result1, result2, result3;
            pystencils::runtime::random::detail::philox_float4(
                ctr0, ctr1, ctr2, ctr3,            //
                (uint32_t)key0, (uint32_t)key1,    //
                result0, result1, result2, result3 //
            );

            return svcreate4(result0, result1, result2, result3);
        }

        template <typename Idx>
        svfloat32x4_t philox_float32(svint32_t ctr0, svint32_t ctr1, svint32_t ctr2, svint32_t ctr3, Idx key0, Idx key1)
        {
            return philox_float32(
                svreinterpret_u32_s32(ctr0),
                svreinterpret_u32_s32(ctr1),
                svreinterpret_u32_s32(ctr2),
                svreinterpret_u32_s32(ctr3),
                key0,
                key1);
        }

        template <typename Idx>
        svfloat64x2_t philox_float64(svuint64_t ctr0, svuint64_t ctr1, svuint64_t ctr2, svuint64_t ctr3, Idx key0, Idx key1)
        {
            svuint32_t ctr0_u32{svuzp1_u32(svreinterpret_u32_u64(ctr0), svdup_n_u32(0u))};
            svuint32_t ctr1_u32{svuzp1_u32(svreinterpret_u32_u64(ctr1), svdup_n_u32(0u))};
            svuint32_t ctr2_u32{svuzp1_u32(svreinterpret_u32_u64(ctr2), svdup_n_u32(0u))};
            svuint32_t ctr3_u32{svuzp1_u32(svreinterpret_u32_u64(ctr3), svdup_n_u32(0u))};

            svfloat64_t result0, result1, ignore;
            pystencils::runtime::random::detail::philox_double2(
                ctr0_u32, ctr1_u32, ctr2_u32, ctr3_u32, //
                (uint32_t)key0, (uint32_t)key1,         //
                result0, ignore, result1, ignore        //
            );

            return svcreate2(result0, result1);
        }
        
	template <typename Idx>
        svfloat64x2_t philox_float64(svint64_t ctr0, svint64_t ctr1, svint64_t ctr2, svint64_t ctr3, Idx key0, Idx key1)
        {
            return philox_float64(
                svreinterpret_u64_s64(ctr0),
                svreinterpret_u64_s64(ctr1),
                svreinterpret_u64_s64(ctr2),
                svreinterpret_u64_s64(ctr3),
                key0,
                key1);
        }

    }
}
