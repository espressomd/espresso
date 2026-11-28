/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#include "for_each_3d.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <iterator>
#include <memory>
#include <new>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

/**
 * @brief Extract a 3D block from the halo field into a caller-provided buffer.
 * @c out must hold at least @c product(stop - start) elements; the block is
 * written in full, so its prior contents are irrelevant.
 */
template <Utils::MemoryOrder memory_order,
          Utils::MemoryOrder output_memory_order, typename Container,
          typename OutValue>
void extract_block_into(OutValue *out, Container const &in_array,
                        Utils::Vector3i const &dimensions,
                        Utils::Vector3i const &start,
                        Utils::Vector3i const &stop) {
  auto const block_dim = stop - start;

  // Extract the block
  if constexpr (memory_order == Utils::MemoryOrder::ROW_MAJOR and
                output_memory_order == Utils::MemoryOrder::ROW_MAJOR) {
    auto const plane_src = dimensions[2] * dimensions[1];
    auto const lane_src = dimensions[2];
    auto const lane_dst = block_dim[2];
    auto const *const src = in_array.data();
    auto *const dst = out;
    auto const n_y = stop[1] - start[1];
    // Explicit destination offset per (x, y) lane lets the outer loop run in
    // parallel: each x writes a disjoint set of contiguous lanes, so the copy
    // is thread-safe and bitwise-identical to the serial version.
    _Pragma("omp parallel for") for (int x = start[0]; x < stop[0]; ++x) {
      for (int y = start[1]; y < stop[1]; ++y) {
        auto const offset_src = x * plane_src + y * lane_src + start[2];
        auto const offset_dst =
            (static_cast<std::size_t>(x - start[0]) * n_y + (y - start[1])) *
            lane_dst;
        std::copy_n(src + offset_src, lane_dst, dst + offset_dst);
      }
    }
  } else {
    for_each_3d_lin<output_memory_order>(
        start, stop, [&](Utils::Vector3i const &indices, int out_index) {
          // Compute indices for input and output arrays
          auto const in_index =
              Utils::get_linear_index<memory_order>(indices, dimensions);
          assert(out_index == Utils::get_linear_index<output_memory_order>(
                                  indices - start, block_dim));
          // Copy the value
          out[out_index] = in_array[in_index];
        });
  }
}

/** @brief Pad a 3D matrix with zeros to restore halo regions, writing into a
 *  caller-provided buffer of @c product(cropped_dim+pad_left+pad_right)
 *  elements. Every element is written exactly once (halo shells are zeroed,
 *  the interior is copied), so the buffer's prior contents are irrelevant and
 *  a persistent buffer can be reused across calls without an allocation or a
 *  full-volume zero-fill.
 */
template <Utils::MemoryOrder memory_order,
          Utils::MemoryOrder output_memory_order, typename T, typename OutValue>
void pad_with_zeros_discard_imag_into(OutValue *out, std::span<T> cropped_array,
                                      Utils::Vector3i const &cropped_dim,
                                      Utils::Vector3i const &pad_left,
                                      Utils::Vector3i const &pad_right) {

  auto constexpr get_real = [](auto const &v) { return std::real(v); };
  // Calculate dimensions and strides
  auto const padded_dim = cropped_dim + pad_left + pad_right;

  if constexpr (memory_order == Utils::MemoryOrder::ROW_MAJOR and
                output_memory_order == Utils::MemoryOrder::ROW_MAJOR) {
    auto const plane_dst = padded_dim[2] * padded_dim[1];
    auto const lane_dst = padded_dim[2];
    auto const lane_src = cropped_dim[2];
    auto *const dst = out;
    auto const *const src = cropped_array.data();
    // One parallel sweep over the padded x-planes: each x writes a disjoint
    // plane (fully zeroed in the x-halo, per-lane zero/copy/zero otherwise),
    // so the sweep is thread-safe and bitwise-identical to the serial version.
    _Pragma("omp parallel for") for (int x = 0; x < padded_dim[0]; ++x) {
      auto *const plane = dst + static_cast<std::size_t>(x) * plane_dst;
      if (x < pad_left[0] or x >= pad_left[0] + cropped_dim[0]) {
        std::fill_n(plane, plane_dst, OutValue{});
        continue;
      }
      auto const src_x = x - pad_left[0];
      for (int y = 0; y < padded_dim[1]; ++y) {
        auto *const lane = plane + static_cast<std::size_t>(y) * lane_dst;
        if (y < pad_left[1] or y >= pad_left[1] + cropped_dim[1]) {
          std::fill_n(lane, lane_dst, OutValue{});
          continue;
        }
        auto const src_y = y - pad_left[1];
        auto const offset_src =
            (static_cast<std::size_t>(src_x) * cropped_dim[1] + src_y) *
            lane_src;
        std::fill_n(lane, pad_left[2], OutValue{});
        if constexpr (std::is_floating_point_v<T>) {
          std::copy_n(src + offset_src, lane_src, lane + pad_left[2]);
        } else {
          std::transform(src + offset_src, src + offset_src + lane_src,
                         lane + pad_left[2], get_real);
        }
        std::fill_n(lane + pad_left[2] + lane_src, pad_right[2], OutValue{});
      }
    }
  } else {
    std::fill_n(out, Utils::product(padded_dim), OutValue{});
    for_each_3d_lin<memory_order>(
        Utils::Vector3i{0, 0, 0}, cropped_dim,
        [&](Utils::Vector3i const &indices, int in_index) {
          // Compute indices for input and output arrays
          auto const out_index = Utils::get_linear_index<output_memory_order>(
              indices + pad_left, padded_dim);
          // Copy the value
          out[out_index] = std::real(cropped_array[in_index]);
        });
  }
}
