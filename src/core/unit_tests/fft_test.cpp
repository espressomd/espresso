/*
 * Copyright (C) 2024 The ESPResSo project
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

#define BOOST_TEST_MODULE "FFT utility functions"
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "config/config.hpp"

#include "fft/fft.hpp"
#include "fft/vector.hpp"
#include "p3m/for_each_3d.hpp"
#include "p3m/packing.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <new>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#if defined(P3M) or defined(DP3M)

BOOST_AUTO_TEST_CASE(fft_find_comm_groups_mismatch) {
  using fft::find_comm_groups;
  int my_pos[3] = {0};
  int nodelist[4] = {0};
  int nodepos[12] = {0};
  int rank = 0;
  {
    auto const optional = find_comm_groups({0, 1, 2}, {1, 2, 3}, nodelist,
                                           nodelist, nodepos, my_pos, rank);
    BOOST_CHECK(not optional.has_value());
  }
  {
    auto const optional = find_comm_groups({3, 2, 1}, {2, 3, 1}, nodelist,
                                           nodelist, nodepos, my_pos, rank);
    BOOST_CHECK(not optional.has_value());
  }
  {
    auto const optional = find_comm_groups({2, 3, 1}, {3, 2, 1}, nodelist,
                                           nodelist, nodepos, my_pos, rank);
    BOOST_CHECK(not optional.has_value());
  }
}

BOOST_AUTO_TEST_CASE(fft_map_grid) {
  using fft::map_3don2d_grid;
  {
    auto g3d = Utils::Vector3i{{3, 2, 1}};
    auto g2d = Utils::Vector3i{{3, 2, 1}};
    auto ref = g2d;
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, 2);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{2, 1, 6}};
    auto g2d = Utils::Vector3i{{6, 2, 1}};
    auto ref = g2d;
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, 2);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{2, 6, 2}};
    auto g2d = Utils::Vector3i{{6, 2, 6}};
    auto ref = Utils::Vector3i{{6, 1, 2}};
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, 1);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{3, 6, 6}};
    auto g2d = Utils::Vector3i{{6, 3, 6}};
    auto ref = g2d;
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, -1);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{4, 1, 6}};
    auto g2d = Utils::Vector3i{{6, 4, 1}};
    auto ref = Utils::Vector3i{{4, 6, 1}};
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, 2);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{5, 7, 7}};
    auto g2d = Utils::Vector3i{{7, 7, 5}};
    auto ref = Utils::Vector3i{{1, 7, 7}};
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, 0);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{5, 7, 5}};
    auto g2d = Utils::Vector3i{{7, 7, 5}};
    auto ref = g2d;
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, -1);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{4, 5, 6}};
    auto g2d = Utils::Vector3i{{6, 4, 5}};
    auto ref = Utils::Vector3i{{4, 1, 6}};
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, 1);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{5, 4, 6}};
    auto g2d = Utils::Vector3i{{6, 4, 5}};
    auto ref = Utils::Vector3i{{1, 4, 6}};
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, 0);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{5, 6, 8}};
    auto g2d = Utils::Vector3i{{8, 7, 5}};
    auto ref = g2d;
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, -1);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
  {
    auto g3d = Utils::Vector3i{{5, 6, 9}};
    auto g2d = Utils::Vector3i{{8, 7, 5}};
    auto ref = g2d;
    auto dir = map_3don2d_grid(g3d.data(), g2d.data());
    BOOST_CHECK_EQUAL(dir, -1);
    BOOST_CHECK_EQUAL(g2d, ref);
  }
}

BOOST_AUTO_TEST_CASE(fft_exceptions) {
  auto constexpr size_max = std::numeric_limits<std::size_t>::max();
  auto constexpr bad_size = size_max / sizeof(float) + 1ul;
  fft::allocator<float> allocator{};
  BOOST_CHECK_EQUAL(allocator.allocate(0ul), nullptr);
  BOOST_CHECK_THROW(allocator.allocate(bad_size), std::bad_array_new_length);
}

BOOST_AUTO_TEST_CASE(fft_plan_without_mpi) {
  auto *plan = new fft::fft_plan<float>();
  plan->destroy_plan(); // no-op since handle is null
  delete plan;
}

#endif // defined(P3M) or defined(DP3M)

BOOST_AUTO_TEST_CASE(for_each_3d_test) {
  auto const m_start = Utils::Vector3i{{0, -1, 3}};
  auto const m_stop = Utils::Vector3i{{2, 2, 5}};
  auto ref_loop_counters = m_start;
  auto indices = Utils::Vector3i{};

  auto const kernel = [&]() {
    BOOST_REQUIRE_EQUAL(indices, ref_loop_counters);
    if (++ref_loop_counters[2u] == m_stop[2u]) {
      ref_loop_counters[2u] = m_start[2u];
      if (++ref_loop_counters[1u] == m_stop[1u]) {
        ref_loop_counters[1u] = m_start[1u];
        if (++ref_loop_counters[0u] == m_stop[0u]) {
          ref_loop_counters[0u] = m_start[0u];
        }
      }
    }
  };

  {
    for_each_3d(m_start, m_stop, indices, kernel, [&](unsigned dim, int n) {
      BOOST_REQUIRE_GE(dim, 0);
      BOOST_REQUIRE_LE(dim, 2);
      BOOST_REQUIRE_EQUAL(n, ref_loop_counters[dim]);
    });

    BOOST_REQUIRE_EQUAL(indices, m_stop);
    BOOST_REQUIRE_EQUAL(ref_loop_counters, m_start);
  }
  {
    for_each_3d(m_start, m_stop, indices, kernel);

    BOOST_REQUIRE_EQUAL(indices, m_stop);
    BOOST_REQUIRE_EQUAL(ref_loop_counters, m_start);
  }
}

template <typename T, int dim1, int dim2, int dim3, typename F>
void run_on_grid(T (&grid)[dim1][dim2][dim3], F fun) {
  for (auto i = 0; i < dim1; ++i) {
    for (auto j = 0; j < dim2; ++j) {
      for (auto k = 0; k < dim3; ++k) {
        fun(grid, i, j, k);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(packing_unpacking) {
  // packing into block
  {
    int const dim[3] = {4, 6, 8};
    int mesh[4][6][8] = {};
    run_on_grid(mesh, [](auto &mesh, int i, int j, int k) {
      mesh[i][j][k] = (i + 1) * (j + 1) * (k + 1);
    });
    int const start[3] = {1, 2, 3};
    int const size[3] = {2, 3, 4};
    int block[2][3][4] = {};
    run_on_grid(block, [](auto &gd, int i, int j, int k) { gd[i][j][k] = -1; });
    fft_pack_block(&(mesh[0][0][0]), &(block[0][0][0]), start, size, dim, 1);
    run_on_grid(block, [&start](auto &block, int i, int j, int k) {
      auto const ref =
          (i + 1 + start[0u]) * (j + 1 + start[1u]) * (k + 1 + start[2u]);
      BOOST_REQUIRE_EQUAL(block[i][j][k], ref);
    });
  }

  // packing into self
  {
    int const dim[3] = {1, 1, 8};
    int mesh[1][1][8] = {};
    run_on_grid(mesh, [](auto &gd, int i, int j, int k) { gd[i][j][k] = k; });
    int const start[3] = {0, 0, 0};
    int const size[3] = {1, 1, 4};
    fft_pack_block(&(mesh[0][0][0]), &(mesh[0][0][3]), start, size, dim, 1);
    run_on_grid(mesh, [&size](auto &mesh, int i, int j, int k) {
      auto ref = k;
      if (k >= 3 and k < 3 + size[2u]) {
        ref = k - 3;
      }
      BOOST_REQUIRE_EQUAL(mesh[i][j][k], ref);
    });
  }

  // unpacking into grid
  {
    int const dim[3] = {4, 6, 8};
    int mesh[4][6][8] = {};
    run_on_grid(mesh, [](auto &mesh, int i, int j, int k) {
      mesh[i][j][k] = (i + 1) * (j + 1) * (k + 1);
    });
    int const start[3] = {1, 2, 3};
    int const size[3] = {2, 3, 4};
    int block[2][3][4] = {};
    run_on_grid(block, [](auto &block, int i, int j, int k) {
      block[i][j][k] = -(i + 1) * (j + 1) * (k + 1);
    });
    fft_unpack_block(&(block[0][0][0]), &(mesh[0][0][0]), start, size, dim, 1);
    run_on_grid(mesh, [&start, &size](auto &mesh, int i, int j, int k) {
      auto ref = (i + 1) * (j + 1) * (k + 1);
      if ((i >= start[0] and i < (start[0] + size[0])) and
          (j >= start[1] and j < (start[1] + size[1])) and
          (k >= start[2] and k < (start[2] + size[2]))) {
        ref = -(i + 1 - start[0]) * (j + 1 - start[1]) * (k + 1 - start[2]);
      }
      BOOST_REQUIRE_EQUAL(mesh[i][j][k], ref);
    });
  }

  // unpacking into self
  {
    int const dim[3] = {1, 1, 8};
    int mesh[1][1][8] = {};
    run_on_grid(mesh, [](auto &gd, int i, int j, int k) { gd[i][j][k] = k; });
    int const start[3] = {0, 0, 0};
    int const size[3] = {1, 1, 4};
    fft_unpack_block(&(mesh[0][0][3]), &(mesh[0][0][0]), start, size, dim, 1);
    run_on_grid(mesh, [&size](auto &mesh, int i, int j, int k) {
      auto const ref = (k < size[2u]) ? k + 3 : k;
      BOOST_REQUIRE_EQUAL(mesh[i][j][k], ref);
    });
  }
}
