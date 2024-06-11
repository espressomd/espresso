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

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

std::optional<std::vector<int>>
find_comm_groups(Utils::Vector3i const &, Utils::Vector3i const &,
                 Utils::Span<const int>, Utils::Span<int>, Utils::Span<int>,
                 Utils::Span<int>, int);

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

#else  // defined(P3M) || defined(DP3M)
int main(int argc, char **argv) {}
#endif // defined(P3M) || defined(DP3M)
