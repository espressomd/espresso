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

#include "config/config.hpp"

#if defined(P3M) || defined(DP3M)

#include <boost/test/unit_test.hpp>

#include "p3m/fft.hpp"
#include "p3m/for_each_3d.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

std::optional<std::vector<int>> find_comm_groups(Utils::Vector3i const &,
                                                 Utils::Vector3i const &,
                                                 std::span<int const>,
                                                 std::span<int>, std::span<int>,
                                                 std::span<int>, int);

BOOST_AUTO_TEST_CASE(fft_find_comm_groups_mismatch) {
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
  auto constexpr bad_size = size_max / sizeof(int) + 1ul;
  fft_allocator<int> allocator{};
  BOOST_CHECK_EQUAL(allocator.allocate(0ul), nullptr);
  BOOST_CHECK_THROW(allocator.allocate(bad_size), std::bad_array_new_length);
}

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

#else  // defined(P3M) || defined(DP3M)
int main(int argc, char **argv) {}
#endif // defined(P3M) || defined(DP3M)
