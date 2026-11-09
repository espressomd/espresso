/*
 * Copyright (C) 2026 The ESPResSo project
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

#define BOOST_TEST_MODULE VoxelMask test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <shapes/VoxelMask.hpp>

#include <utils/Vector.hpp>

#include <cmath>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <vector>

using Shapes::VoxelMask;

namespace {
VoxelMask::Geometry default_geometry() {
  VoxelMask::Geometry g;
  g.mask_shape = {10, 10, 10};
  g.voxel_size = 1.;
  g.origin = {0., 0., 0.};
  g.box_l = {10., 10., 10.};
  g.periodic = {true, true, true};
  g.halo = 2;
  return g;
}

Utils::Vector3d const local_lower{0., 0., 0.};
Utils::Vector3d const local_upper{5., 5., 5.};

/** Fill the stored region from a predicate on mask indices. */
void load(VoxelMask &mask, std::function<bool(Utils::Vector3i const &)> solid) {
  auto const &g = mask.geometry();
  auto const region = VoxelMask::stored_region(g, local_lower, local_upper);
  std::vector<std::uint8_t> voxels;
  for (int i = 0; i < region.extent[0]; ++i) {
    for (int j = 0; j < region.extent[1]; ++j) {
      for (int k = 0; k < region.extent[2]; ++k) {
        Utils::Vector3i const unfolded{region.lower[0] + i, region.lower[1] + j,
                                       region.lower[2] + k};
        Utils::Vector3i m;
        bool pore = false;
        for (int a = 0; a < 3; ++a) {
          m[a] = VoxelMask::mask_index(g, a, unfolded[a]);
          pore = pore or m[a] < 0;
        }
        voxels.push_back((not pore and solid(m)) ? 1 : 0);
      }
    }
  }
  mask.set_local_data(region, voxels);
}

void check_close(Utils::Vector3d const &a, Utils::Vector3d const &b) {
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_SMALL(a[i] - b[i], 1e-12);
  }
}
} // namespace

BOOST_AUTO_TEST_CASE(stored_region_bounds) {
  auto const region =
      VoxelMask::stored_region(default_geometry(), local_lower, local_upper);
  for (int a = 0; a < 3; ++a) {
    BOOST_CHECK_EQUAL(region.lower[a], -2);
    BOOST_CHECK_EQUAL(region.extent[a], 9);
  }
}

BOOST_AUTO_TEST_CASE(mask_index_periodic_wrap) {
  auto g = default_geometry();
  BOOST_CHECK_EQUAL(VoxelMask::mask_index(g, 0, -1), 9);
  BOOST_CHECK_EQUAL(VoxelMask::mask_index(g, 0, 10), 0);
  g.periodic[0] = false;
  BOOST_CHECK_EQUAL(VoxelMask::mask_index(g, 0, -1), -1);
  BOOST_CHECK_EQUAL(VoxelMask::mask_index(g, 0, 10), -1);
}

BOOST_AUTO_TEST_CASE(mask_index_box_not_multiple_of_voxel) {
  auto g = default_geometry();
  g.box_l = {10.5, 10.5, 10.5};
  g.mask_shape = {11, 11, 11};
  BOOST_CHECK_EQUAL(VoxelMask::mask_index(g, 0, 9), 9);
  BOOST_CHECK_EQUAL(VoxelMask::mask_index(g, 0, 10), 0);
  BOOST_CHECK_EQUAL(VoxelMask::mask_index(g, 0, 11), 1);
}

BOOST_AUTO_TEST_CASE(slab_distance) {
  VoxelMask mask;
  mask.set_geometry(default_geometry());
  load(mask, [](auto const &m) { return m[0] < 2; });
  double dist;
  Utils::Vector3d vec;

  Utils::Vector3d pos{3.5, 2., 2.};
  mask.calculate_dist(pos, dist, vec);
  BOOST_CHECK_CLOSE(dist, 1.5, 1e-9);
  check_close(vec, {1.5, 0., 0.});
  BOOST_CHECK(not mask.is_inside(pos));

  pos = {1.5, 2., 2.};
  mask.calculate_dist(pos, dist, vec);
  BOOST_CHECK_CLOSE(dist, -0.5, 1e-9);
  check_close(vec, {-0.5, 0., 0.});
  BOOST_CHECK(mask.is_inside(pos));
}

BOOST_AUTO_TEST_CASE(edge_distance) {
  VoxelMask mask;
  mask.set_geometry(default_geometry());
  load(mask, [](auto const &m) { return m == Utils::Vector3i{2, 2, 2}; });
  double dist;
  Utils::Vector3d vec;
  mask.calculate_dist({3.5, 3.5, 2.5}, dist, vec);
  BOOST_CHECK_CLOSE(dist, std::sqrt(0.5), 1e-9);
  check_close(vec, {0.5, 0.5, 0.});
}

BOOST_AUTO_TEST_CASE(far_from_solid) {
  VoxelMask mask;
  mask.set_geometry(default_geometry());
  load(mask, [](auto const &) { return false; });
  double dist;
  Utils::Vector3d vec;
  mask.calculate_dist({2.5, 2.5, 2.5}, dist, vec);
  // beyond any interaction cutoff, with a finite square
  BOOST_CHECK_GE(dist, 1e100);
  BOOST_CHECK(std::isfinite(dist * dist));
  check_close(vec / dist, {1., 0., 0.});
  BOOST_CHECK_EQUAL(mask.local_distance_reach(), 2.);

  load(mask, [](auto const &) { return true; });
  mask.calculate_dist({2.5, 2.5, 2.5}, dist, vec);
  BOOST_CHECK_LE(dist, -1e100);
  BOOST_CHECK(std::isfinite(dist * dist));
  check_close(vec / dist, {1., 0., 0.});
}

BOOST_AUTO_TEST_CASE(image_position_across_periodic_face) {
  VoxelMask mask;
  mask.set_geometry(default_geometry());
  load(mask, [](auto const &m) { return m[0] == 9; });
  double dist;
  Utils::Vector3d vec;

  Utils::Vector3d pos{9.8, 2.5, 2.5};
  BOOST_CHECK(mask.is_inside(pos));
  mask.calculate_dist(pos, dist, vec);
  BOOST_CHECK_CLOSE(dist, -0.2, 1e-9);
  check_close(vec, {-0.2, 0., 0.});

  mask.calculate_dist({0.3, 2.5, 2.5}, dist, vec);
  BOOST_CHECK_CLOSE(dist, 0.3, 1e-9);
  check_close(vec, {0.3, 0., 0.});
}

BOOST_AUTO_TEST_CASE(outside_region_throws) {
  VoxelMask mask;
  mask.set_geometry(default_geometry());
  load(mask, [](auto const &) { return false; });
  double dist;
  Utils::Vector3d vec;
  BOOST_CHECK_THROW(mask.calculate_dist({7.5, 2., 2.}, dist, vec),
                    std::runtime_error);
  BOOST_CHECK_THROW(mask.is_inside({7.5, 2., 2.}), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(memory_is_bit_packed) {
  VoxelMask mask;
  mask.set_geometry(default_geometry());
  load(mask, [](auto const &m) { return m[0] < 2; });
  BOOST_CHECK_EQUAL(mask.local_voxel_count(), 9u * 9u * 9u);
  BOOST_CHECK_LE(mask.local_memory_bytes(), mask.local_voxel_count() / 8 + 64);
}

BOOST_AUTO_TEST_CASE(loader_called_on_geometry_change) {
  VoxelMask mask;
  int calls = 0;
  mask.set_loader([&calls](VoxelMask &) { ++calls; });
  mask.on_local_geometry_change();
  BOOST_CHECK_EQUAL(calls, 1);
  mask.on_local_geometry_change();
  BOOST_CHECK_EQUAL(calls, 2);
}

BOOST_AUTO_TEST_CASE(block_wise_fill_matches_whole_region) {
  auto const g = default_geometry();
  auto const region = VoxelMask::stored_region(g, local_lower, local_upper);
  auto const solid = [](int i, int j, int k) { return (i + 2 * j + k) % 3; };
  auto const &e = region.extent;
  std::vector<std::uint8_t> all;
  for (int i = 0; i < e[0]; ++i)
    for (int j = 0; j < e[1]; ++j)
      for (int k = 0; k < e[2]; ++k)
        all.push_back(solid(i, j, k) == 0 ? 1 : 0);
  VoxelMask whole;
  whole.set_geometry(g);
  whole.set_local_data(region, all);

  // fill by x slabs of two planes, the last slab split into rows
  VoxelMask blocks;
  blocks.set_geometry(g);
  blocks.allocate_local_data(region);
  BOOST_CHECK_EQUAL(blocks.local_voxel_count(), whole.local_voxel_count());
  auto const fill = [&](Utils::Vector3i const &lower,
                        Utils::Vector3i const &extent) {
    std::vector<std::uint8_t> voxels;
    for (int i = 0; i < extent[0]; ++i)
      for (int j = 0; j < extent[1]; ++j)
        for (int k = 0; k < extent[2]; ++k)
          voxels.push_back(
              solid(lower[0] + i, lower[1] + j, lower[2] + k) == 0 ? 1 : 0);
    blocks.set_local_block(lower, extent, voxels);
  };
  int x = 0;
  for (; x + 2 <= e[0]; x += 2) {
    fill({x, 0, 0}, {2, e[1], e[2]});
  }
  for (int y = 0; y < e[1]; ++y) {
    fill({x, y, 0}, {e[0] - x, 1, e[2]});
  }
  for (int i = 0; i < e[0]; ++i)
    for (int j = 0; j < e[1]; ++j)
      for (int k = 0; k < e[2]; ++k) {
        Utils::Vector3d const pos{i - 1.5, j - 1.5, k - 1.5};
        BOOST_REQUIRE_EQUAL(blocks.is_inside(pos), whole.is_inside(pos));
      }

  std::vector<std::uint8_t> const too_small(3);
  BOOST_CHECK_THROW(blocks.set_local_block({0, 0, 0}, {1, 1, 4}, too_small),
                    std::runtime_error);
  std::vector<std::uint8_t> const one(1);
  BOOST_CHECK_THROW(blocks.set_local_block({e[0], 0, 0}, {1, 1, 1}, one),
                    std::runtime_error);
  BOOST_CHECK_THROW(blocks.set_local_block({-1, 0, 0}, {1, 1, 1}, one),
                    std::runtime_error);
}
