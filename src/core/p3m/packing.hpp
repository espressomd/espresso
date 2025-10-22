/*
 * Copyright (C) 2010-2024 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include <cstddef>
#include <cstring>

/**
 * @brief Pack a 3D-block of size @c size starting at @c start of an input
 * 3D-grid @c in with dimension @c dim into an output 3D-block @c out.
 *
 * The block is stored in 'row-major-order' or 'C-order', that means
 * the first index is changing slowest when running through the linear array.
 * The element (@c i0 (slow), @c i1 (mid), @c i2 (fast)) has the linear
 * index <tt>li = i2 + size[2] * (i1 + (size[1] * i0))</tt>.
 *
 * @param[in]  in      pointer to input 3D-grid.
 * @param[out] out     pointer to output 3D-block.
 * @param[in]  start   start index of the block in the input 3D-grid.
 * @param[in]  size    size of the output 3D-block.
 * @param[in]  dim     size of the input 3D-grid.
 * @param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
template <typename FloatType>
void fft_pack_block(FloatType const *const in, FloatType *const out,
                    int const *start, int const *size, int const *dim,
                    int element) {

  auto const copy_size =
      static_cast<std::size_t>(element * size[2]) * sizeof(FloatType);
  /* offsets for indices in input grid */
  auto const m_in_offset = element * dim[2];
  auto const s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  /* offsets for indices in output grid */
  auto const m_out_offset = element * size[2];
  /* linear index of input grid, linear index of output grid */
  int li_in = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));
  int li_out = 0;

  for (int s = 0; s < size[0]; s++) {
    for (int m = 0; m < size[1]; m++) {
      memmove(&(out[li_out]), &(in[li_in]), copy_size);
      li_in += m_in_offset;
      li_out += m_out_offset;
    }
    li_in += s_in_offset;
  }
}

/**
 * @brief Unpack a 3D-block @c in of size @c size into an output
 * 3D-grid @c out of size @c dim starting at position @c start.
 *
 * See also @ref fft_pack_block.
 *
 * @param[in]  in      pointer to input 3D-block.
 * @param[out] out     pointer to output 3D-grid.
 * @param[in]  start   start index of the block in the input 3D-grid.
 * @param[in]  size    size of the input 3D-block.
 * @param[in]  dim     size of the output 3D-grid.
 * @param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
template <typename FloatType>
void fft_unpack_block(FloatType const *const in, FloatType *const out,
                      int const *start, int const *size, int const *dim,
                      int element) {

  auto const copy_size =
      static_cast<std::size_t>(element * size[2]) * sizeof(FloatType);
  /* offsets for indices in output grid */
  auto const m_out_offset = element * dim[2];
  auto const s_out_offset = element * (dim[2] * (dim[1] - size[1]));
  /* offset for indices in input grid */
  auto const m_in_offset = element * size[2];
  /* linear index of in grid, linear index of out grid */
  int li_in = 0;
  int li_out = element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (int s = 0; s < size[0]; s++) {
    for (int m = 0; m < size[1]; m++) {
      memmove(&(out[li_out]), &(in[li_in]), copy_size);
      li_in += m_in_offset;
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}
