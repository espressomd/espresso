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

#include <cassert>
#include <cstddef>
#include <span>
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

  return out_array;
}

// Function to pad the 3D cropped field with zeros to restore the halo regions
template <Utils::MemoryOrder memory_order,
          Utils::MemoryOrder output_memory_order, typename T>
auto pad_with_zeros_discard_imag(std::span<T> cropped_array,
                                 Utils::Vector3i const &cropped_dim,
                                 Utils::Vector3i const &pad_left,
                                 Utils::Vector3i const &pad_right) {

  // Calculate dimensions and strides
  Utils::Vector3i padded_dim = cropped_dim + pad_left + pad_right;
  // Output vector to hold the padded field (initialized with zeros)
  std::vector<typename T::value_type> padded_array(Utils::product(padded_dim));

  // Fill in the original cropped field into the padded array by chunks
  for_each_3d_lin<memory_order>(
      Utils::Vector3i{0, 0, 0}, cropped_dim,
      [&](Utils::Vector3i const &indices, int in_index) {
        // Compute indices for input and output arrays
        auto const out_index = Utils::get_linear_index<output_memory_order>(
            indices + pad_left, padded_dim);
        // Copy the value
        padded_array[out_index] = cropped_array[in_index].real();
      });

  return padded_array;
}
