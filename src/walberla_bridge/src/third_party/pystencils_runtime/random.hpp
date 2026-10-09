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
#include "./short_arrays.hpp"
#include "./bits/philox_rand.h"

#if PYSTENCILS_PLATFORM_IS_HIP
#include "./hip.h"
#endif

namespace pystencils::runtime::random
{

    template <typename Idx>
    PYSTENCILS_DEVICE_QUALIFIERS ShortArray<float, 4> philox_float32(Idx ctr0, Idx ctr1, Idx ctr2, Idx ctr3, Idx key0, Idx key1)
    {
        ShortArray<float, 4> result;
        detail::philox_float4(
            (uint32_t)ctr0, (uint32_t)ctr1, (uint32_t)ctr2, (uint32_t)ctr3, //
            (uint32_t)key0, (uint32_t)key1,                                 //
            result[0], result[1], result[2], result[3]                      //
        );
        return result;
    }

    template <typename Idx>
    PYSTENCILS_DEVICE_QUALIFIERS ShortArray<double, 2> philox_float64(Idx ctr0, Idx ctr1, Idx ctr2, Idx ctr3, Idx key0, Idx key1)
    {
        ShortArray<double, 2> result;
        detail::philox_double2(
            (uint32_t)ctr0, (uint32_t)ctr1, (uint32_t)ctr2, (uint32_t)ctr3, //
            (uint32_t)key0, (uint32_t)key1,                                 //
            result[0], result[1]                                            //
        );
        return result;
    }

} // pystencils::runtime::random

#undef QUALIFIERS
