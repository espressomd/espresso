/*
 * Copyright (C) 2024-2025 The ESPResSo project
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
#include <span>
#include <type_traits>
#include <vector>

// Function to extract a 3D block from the halo field
template <Utils::MemoryOrder memory_order,
          Utils::MemoryOrder output_memory_order, typename Container>
auto extract_block(Container const &in_array, Utils::Vector3i const &dimensions,
                   Utils::Vector3i const &start, Utils::Vector3i const &stop) {
  // Calculate the size of the block excluding halo regions
  auto const block_dim = stop - start;
  auto const size = static_cast<std::size_t>(Utils::product(block_dim));

  // Output vector to hold the block
  std::vector<typename Container::value_type> out_array(size);

  // Extract the block
  if constexpr (memory_order == Utils::MemoryOrder::ROW_MAJOR and
                output_memory_order == Utils::MemoryOrder::ROW_MAJOR) {
    auto const plane_src = dimensions[2] * dimensions[1];
    auto const lane_src = dimensions[2];
    auto const lane_dst = block_dim[2];
    auto src = in_array.data();
    auto dst = out_array.begin();
    for (int x = start[0]; x < stop[0]; ++x) {
      for (int y = start[1]; y < stop[1]; ++y) {
        auto const offset_src = x * plane_src + y * lane_src + start[2];
        std::copy_n(src + offset_src, lane_dst, dst);
        std::advance(dst, lane_dst);
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
          out_array[out_index] = in_array[in_index];
        });
  }

  return out_array;
}

/** @brief Pad a 3D matrix with zeros to restore halo regions. */
template <Utils::MemoryOrder memory_order,
          Utils::MemoryOrder output_memory_order, typename T>
auto pad_with_zeros_discard_imag(std::span<T> cropped_array,
                                 Utils::Vector3i const &cropped_dim,
                                 Utils::Vector3i const &pad_left,
                                 Utils::Vector3i const &pad_right) {

  using value_type = decltype(std::real(std::declval<T>()));
  auto constexpr get_real = [](auto const &v) { return std::real(v); };
  // Calculate dimensions and strides
  auto const padded_dim = cropped_dim + pad_left + pad_right;
  // Output vector to hold the padded field (initialized with zeros)
  std::vector<value_type> padded_array(Utils::product(padded_dim));

  // Fill in the original cropped field into the padded array by chunks
  if constexpr (memory_order == Utils::MemoryOrder::ROW_MAJOR and
                output_memory_order == Utils::MemoryOrder::ROW_MAJOR) {
    auto const plane_dst = padded_dim[2] * padded_dim[1];
    auto const lane_dst = padded_dim[2];
    auto const lane_src = cropped_dim[2];
    auto const base_offset_dst =
        pad_left[0] * plane_dst + pad_left[1] * lane_dst + pad_left[2];
    auto dst = padded_array.data();
    auto src = cropped_array.begin();
    for (int x = 0; x < cropped_dim[0]; ++x) {
      for (int y = 0; y < cropped_dim[1]; ++y) {
        auto const offset_dst = base_offset_dst + x * plane_dst + y * lane_dst;
        if constexpr (std::is_floating_point_v<T>) {
          std::copy_n(src, lane_src, dst + offset_dst);
        } else {
          std::transform(src, src + lane_src, dst + offset_dst, get_real);
        }
        std::advance(src, lane_src);
      }
    }
  } else {
    for_each_3d_lin<memory_order>(
        Utils::Vector3i{0, 0, 0}, cropped_dim,
        [&](Utils::Vector3i const &indices, int in_index) {
          // Compute indices for input and output arrays
          auto const out_index = Utils::get_linear_index<output_memory_order>(
              indices + pad_left, padded_dim);
          // Copy the value
          padded_array[out_index] = std::real(cropped_array[in_index]);
        });
  }

  return padded_array;
}
