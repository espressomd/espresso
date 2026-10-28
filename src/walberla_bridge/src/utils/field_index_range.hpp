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
 * Check of the size of the fields of a block, for the kernels that compute
 * their offsets in the field with 32 bits.
 */

#include <core/DataTypes.h>
#include <field/allocation/FieldAllocator.h>

#include <cstddef>
#include <stdexcept>
#include <string>

namespace walberla::detail {

/**
 * @brief Check that a field of a block can be indexed with 32 bits.
 *
 * Several generated kernels compute all their offsets in a field with
 * @c int32_t, for every precision and target: the link kernels of the
 * boundary conditions (@c DynamicUBB of the LB, @c Dirichlet and @c FixedFlux
 * of the EK) and the single-precision AVX stream-collide kernels
 * (@c index_dtype=int32). The number of elements of the field of a block,
 * ghost layers included, must therefore stay below 2^31, otherwise the
 * offsets overflow silently. The kernels don't check this, so it is checked
 * by the constructor of the class that allocates the fields, before it
 * allocates any field.
 *
 * With the aligned allocator (@c field::AllocateAligned, used by the AVX
 * kernels), waLBerla pads the x extent of the field with ghost layers (the
 * fastest coordinate in the @c fzyx layout) up to a multiple of
 * <tt>alignment / sizeof(T)</tt> elements (see @c allocateMemory in
 * <tt>src/field/allocation/FieldAllocator.h</tt>), and the strides of the
 * kernels derive from the padded allocation. The padding is applied here
 * unconditionally, to the largest multiple the bridge can use: 16 elements,
 * i.e. 64 bytes (AVX-512) for the smallest element type (single precision).
 * This overestimates the size by less than 1% for any block that can
 * approach the limit, and keeps the check independent of the build options.
 *
 * @param field_name Description of the field for the error message
 *                   (e.g. "PDF field")
 * @param kernels Description of the kernels that index the field for the
 *                error message (e.g. "LB kernels")
 * @param n_components Number of components per cell (e.g. @c Stencil::Size)
 * @param cells_x Number of cells of a block in x direction (without halo)
 * @param cells_y Number of cells of a block in y direction (without halo)
 * @param cells_z Number of cells of a block in z direction (without halo)
 * @param n_ghost_layers Number of ghost layers on each side of a block
 * @throw std::runtime_error if the (padded) field has 2^31 elements or more
 */
inline void check_field_index_range(std::string const &field_name,
                                    std::string const &kernels,
                                    std::size_t n_components, uint_t cells_x,
                                    uint_t cells_y, uint_t cells_z,
                                    uint_t n_ghost_layers) {
  // the factors are checked one at a time, so that nothing can overflow
  auto constexpr limit = uint_t{1} << 31u;
  // padding of the x extent: 64 bytes in elements of 4 bytes
  auto constexpr padding = uint_t{16};
  static_assert(field::SIMDAlignment() <= padding * sizeof(float),
                "The padding of the aligned fields exceeds the assumption");
  auto size = static_cast<uint_t>(n_components);
  auto is_x = true;
  for (auto const cells : {cells_x, cells_y, cells_z}) {
    auto extent = cells + 2u * n_ghost_layers;
    if (is_x and extent <= limit) {
      extent = (extent + padding - 1u) / padding * padding;
    }
    is_x = false;
    if (extent != 0u and size > (limit - 1u) / extent) {
      throw std::runtime_error(
          "The " + field_name + " of a block (" + std::to_string(n_components) +
          " components in " + std::to_string(cells_x) + "x" +
          std::to_string(cells_y) + "x" + std::to_string(cells_z) +
          " cells with " + std::to_string(n_ghost_layers) +
          " ghost layers) has 2^31 or more elements, which the 32-bit index "
          "arithmetic of the " +
          kernels +
          " cannot address. Use more MPI ranks or blocks, or a smaller grid.");
    }
    size *= extent;
  }
}

} // namespace walberla::detail
