/*
 * Copyright (C) 2026 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

/**
 * @file
 * Force-included into the generated CUDA kernels when compiling with Clang.
 *
 * The pystencils runtime marks its device functions @c __device__ only when
 * @c __CUDA_ARCH__ is defined, i.e. during the device compilation pass.
 * Clang also type-checks @c __global__ kernel bodies during the host pass,
 * where these functions are then host-only and cannot be called.
 * Apply the rule of @c pystencils_runtime/bits/philox_rand.h instead,
 * which marks its functions @c __device__ in both Clang passes.
 */

#include <pystencils_runtime/bits/platforms.h>

#if defined(__clang__) && defined(__CUDA__) && !defined(__CUDA_ARCH__)
#undef PYSTENCILS_DEVICE_QUALIFIERS
#define PYSTENCILS_DEVICE_QUALIFIERS __device__ __forceinline__
#endif
