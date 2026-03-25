/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#include <config/config.hpp>

#include "p3m/field_layout_helpers.hpp"

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
#include <Kokkos_Core.hpp>
#endif

#include <utils/Vector.hpp>
#include <utils/index.hpp>

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

struct GlobalConfig {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
#endif
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);
BOOST_AUTO_TEST_SUITE(suite)

template <Utils::MemoryOrder MemOrderReal, Utils::MemoryOrder MemOrderFourier>
void check_add_remove_halo() {
  auto constexpr n_x = 3;
  auto constexpr n_y = 7;
  auto constexpr n_z = 11;
  Utils::Vector3i const shape = {n_x, n_y, n_z};
  Utils::Vector3i const halo_left = {1, 2, 3};
  Utils::Vector3i const halo_right = {3, 1, 2};
  auto const buf_size = static_cast<std::size_t>(Utils::product(shape));
  std::vector<double> real_without_halo_orig(buf_size);
  std::vector<std::complex<double>> cplx_without_halo_orig(buf_size);

  // Fill
  for (int i = 0; i < n_x; ++i) {
    for (int j = 0; j < n_y; ++j) {
      for (int k = 0; k < n_z; ++k) {
        auto const ind = Utils::get_linear_index<MemOrderReal>(i, j, k, shape);
        real_without_halo_orig[ind] = static_cast<double>(ind);
        cplx_without_halo_orig[ind] = {static_cast<double>(ind), 1.};
      }
    }
  }

  // add halo
  auto const shape_with_halo = shape + halo_left + halo_right;
  auto const buf_size_with_halo =
      static_cast<std::size_t>(Utils::product(shape_with_halo));
  auto const real_with_halo =
      pad_with_zeros_discard_imag<MemOrderReal, MemOrderFourier>(
          std::span(real_without_halo_orig), shape, halo_left, halo_right);
  auto const cplx_with_halo =
      pad_with_zeros_discard_imag<MemOrderReal, MemOrderFourier>(
          std::span(cplx_without_halo_orig), shape, halo_left, halo_right);

  BOOST_CHECK_EQUAL(real_with_halo.size(), buf_size_with_halo);
  BOOST_CHECK_EQUAL(cplx_with_halo.size(), buf_size_with_halo);

  // verify inner region
  for (int i = 0; i < n_x; ++i) {
    for (int j = 0; j < n_y; ++j) {
      for (int k = 0; k < n_z; ++k) {
        auto const ind_without_halo =
            Utils::get_linear_index<MemOrderReal>(i, j, k, shape);
        auto const ind_with_halo = Utils::get_linear_index<MemOrderFourier>(
            i + halo_left[0], j + halo_left[1], k + halo_left[2],
            shape_with_halo);
        BOOST_CHECK_EQUAL(real_with_halo[ind_with_halo],
                          real_without_halo_orig[ind_without_halo]);
        BOOST_CHECK_EQUAL(cplx_with_halo[ind_with_halo],
                          cplx_without_halo_orig[ind_without_halo].real());
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
          auto const ind = Utils::get_linear_index<MemOrderFourier>(
              i, j, k, shape_with_halo);
          BOOST_CHECK_EQUAL(real_with_halo[ind], 0.);
          BOOST_CHECK_EQUAL(cplx_with_halo[ind], 0.);
        }
      }
    }
  }

  // Remove halo again
  auto real_without_halo_new = extract_block<MemOrderFourier, MemOrderReal>(
      real_with_halo, shape_with_halo, halo_left, shape_with_halo - halo_right);
  auto cplx_without_halo_new = extract_block<MemOrderFourier, MemOrderReal>(
      cplx_with_halo, shape_with_halo, halo_left, shape_with_halo - halo_right);
  for (std::size_t i = 0u; i < buf_size; ++i) {
    BOOST_CHECK_EQUAL(real_without_halo_new[i], real_without_halo_orig[i]);
    BOOST_CHECK_EQUAL(cplx_without_halo_new[i],
                      cplx_without_halo_orig[i].real());
  }
  BOOST_CHECK_EQUAL(real_without_halo_new.size(), buf_size);
  BOOST_CHECK_EQUAL(cplx_without_halo_new.size(), buf_size);
}

BOOST_AUTO_TEST_CASE(add_remove_halo_col_row) {
  check_add_remove_halo<Utils::MemoryOrder::COLUMN_MAJOR,
                        Utils::MemoryOrder::ROW_MAJOR>();
}
BOOST_AUTO_TEST_CASE(add_remove_halo_col_col) {
  check_add_remove_halo<Utils::MemoryOrder::COLUMN_MAJOR,
                        Utils::MemoryOrder::COLUMN_MAJOR>();
}
BOOST_AUTO_TEST_CASE(add_remove_halo_row_row) {
  check_add_remove_halo<Utils::MemoryOrder::ROW_MAJOR,
                        Utils::MemoryOrder::ROW_MAJOR>();
}
BOOST_AUTO_TEST_CASE(add_remove_halo_row_col) {
  check_add_remove_halo<Utils::MemoryOrder::ROW_MAJOR,
                        Utils::MemoryOrder::COLUMN_MAJOR>();
}

BOOST_AUTO_TEST_SUITE_END()
