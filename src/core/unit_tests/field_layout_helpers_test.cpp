/*
 * Copyright (C) 2024-2025 The ESPResSo project
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

#define BOOST_TEST_MODULE "P3M field layout functions"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "p3m/field_layout_helpers.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

BOOST_AUTO_TEST_CASE(add_remove_halo) {
  int n_x = 3;
  int n_y = 7;
  int n_z = 11;
  std::vector<std::complex<double>> without_halo_orig(n_x * n_y * n_z);
  Utils::Vector3i const shape = {n_x, n_y, n_z};
  Utils::Vector3i const halo_left = {1, 2, 3};
  Utils::Vector3i const halo_right = {3, 1, 2};

  // Fill
  for (int i = 0; i < n_x; ++i) {
    for (int j = 0; j < n_y; ++j) {
      for (int k = 0; k < n_z; ++k) {
        auto const ind = Utils::get_linear_index(i, j, k, shape,
                                                 Utils::MemoryOrder::ROW_MAJOR);
        without_halo_orig[ind] = {static_cast<double>(ind), 0.};
      }
    }
  }

  // add halo
  auto const shape_with_halo = shape + halo_left + halo_right;
  auto const with_halo =
      pad_with_zeros_discard_imag<Utils::MemoryOrder::ROW_MAJOR,
                                  Utils::MemoryOrder::ROW_MAJOR>(
          std::span{without_halo_orig.data(), without_halo_orig.size()}, shape,
          halo_left, halo_right);

  BOOST_CHECK_EQUAL(with_halo.size(), shape_with_halo[0] * shape_with_halo[1] *
                                          shape_with_halo[2]);

  // verify inner region
  for (int i = 0; i < n_x; ++i) {
    for (int j = 0; j < n_y; ++j) {
      for (int k = 0; k < n_z; ++k) {
        auto const ind_without_halo = Utils::get_linear_index(
            i, j, k, shape, Utils::MemoryOrder::ROW_MAJOR);
        auto const ind_with_halo = Utils::get_linear_index(
            i + halo_left[0], j + halo_left[1], k + halo_left[2],
            shape_with_halo, Utils::MemoryOrder::ROW_MAJOR);
        BOOST_CHECK_EQUAL(with_halo[ind_with_halo],
                          without_halo_orig[ind_without_halo]);
      }
    }
  }

  // verify zero padding
  for (int i = 0; i < shape_with_halo[0]; ++i) {
    for (int j = 0; j < shape_with_halo[1]; ++j) {
      for (int k = 0; k < shape_with_halo[2]; ++k) {
        if ((i < halo_left[0] or i >= shape_with_halo[1] - halo_right[0]) or
            (j < halo_left[1] or j >= shape_with_halo[1] - halo_right[1]) or
            (k < halo_left[2] or k >= shape_with_halo[2] - halo_right[2])) {
          // this is a halo cell, must be zero
          auto const ind = Utils::get_linear_index(
              i, j, k, shape_with_halo, Utils::MemoryOrder::ROW_MAJOR);
          BOOST_CHECK_EQUAL(with_halo[ind], 0.);
        }
      }
    }
  }

  // Remove halo again
  auto without_halo_new = extract_block<Utils::MemoryOrder::ROW_MAJOR,
                                        Utils::MemoryOrder::ROW_MAJOR>(
      with_halo, shape_with_halo, halo_left, shape_with_halo - halo_right);
  for (auto i = std::size_t{0u}; i < without_halo_orig.size(); ++i) {
    BOOST_CHECK_EQUAL(without_halo_new[i], without_halo_orig[i]);
  }
  BOOST_CHECK_EQUAL(without_halo_new.size(), without_halo_orig.size());
}
