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

#include "config/config.hpp"

#include "for_each_3d.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <cstddef>
#include <span>
#include <vector>

// Function to extract a 3D block from the halo field
template <Utils::MemoryOrder memory_order_in,
          Utils::MemoryOrder memory_order_out, typename Container>
auto extract_block(Container const &in_array, Utils::Vector3i const &dimensions,
                   Utils::Vector3i const &start, Utils::Vector3i const &stop) {
  // Calculate the size of the block excluding halo regions
  auto const block_dim = stop - start;
  auto const size = static_cast<std::size_t>(Utils::product(block_dim));

  // Output vector to hold the block
  std::vector<typename Container::value_type> out_array(size);

  // Extract the block
  std::size_t out_index{};
  Utils::Vector3i indices{};
  for_each_3d_lin<memory_order_out>(start, stop, indices, out_index, [&]() {
    auto const in_index =
        Utils::get_linear_index<memory_order_in>(indices, dimensions);
#ifdef ADDITIONAL_CHECKS
    assert(out_index == Utils::get_linear_index<memory_order_out>(
                            indices - start, block_dim));
#endif
    out_array[out_index] = in_array[in_index];
  });

  return out_array;
}

// Function to pad the 3D cropped field with zeros to restore the halo regions
template <typename T>
auto pad_with_zeros_discard_imag(std::span<T> cropped_array,
                                 Utils::Vector3i cropped_dim,
                                 Utils::Vector3i pad_left,
                                 Utils::Vector3i pad_right) {

  // Calculate dimensions and strides
  Utils::Vector3i padded_dim = cropped_dim + pad_left + pad_right;
  auto const cropped_xy_stride = cropped_dim[1] * cropped_dim[2];
  auto const padded_xy_stride = padded_dim[1] * padded_dim[2];

  // Output vector to hold the padded field (initialized with zeros)
  std::vector<typename T::value_type> padded_array(
      padded_dim[0] * padded_dim[1] * padded_dim[2]);

  // Calculate the starting position in the padded array for the inner field
  auto const padded_start_x = pad_left[0] * padded_xy_stride;
  auto const padded_start_y = pad_left[1] * padded_dim[2] + pad_left[2];

  // Fill in the original cropped field into the padded array by chunks
  for (int x = 0; x < cropped_dim[0]; ++x) {
    auto const cropped_x_offset = x * cropped_xy_stride;
    auto const padded_x_offset = padded_start_x + x * padded_xy_stride;

    for (int y = 0; y < cropped_dim[1]; ++y) {
      auto const cropped_y_offset = cropped_x_offset + y * cropped_dim[2];
      auto const padded_y_offset =
          padded_x_offset + y * padded_dim[2] + padded_start_y;

      // Copy a contiguous slice of the z-dimension at once
      for (int i = 0; i < cropped_dim[2]; i++) {
        padded_array[padded_y_offset + i] =
            cropped_array[cropped_y_offset + i].real();
      }
    }
  }

  return padded_array;
}
