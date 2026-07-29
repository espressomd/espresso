/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BoxGeometry.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

BOOST_AUTO_TEST_CASE(periodicity_test) {
  /* getter/default */
  {
    auto const box = BoxGeometry{};

    BOOST_CHECK(box.periodic(0u));
    BOOST_CHECK(box.periodic(1u));
    BOOST_CHECK(box.periodic(2u));
  }

  /* setter */
  {
    auto box = BoxGeometry{};

    box.set_periodic(0u, false);
    BOOST_CHECK(not box.periodic(0u));
    box.set_periodic(1u, false);
    BOOST_CHECK(not box.periodic(1u));
    box.set_periodic(2u, false);
    BOOST_CHECK(not box.periodic(2u));
    BOOST_CHECK_THROW(box.set_periodic(3u, false), std::out_of_range);
  }
}

BOOST_AUTO_TEST_CASE(length_test) {
  /* getter/default */
  {
    auto const box = BoxGeometry{};

    BOOST_CHECK(Utils::Vector3d::broadcast(1.) == box.length());
    BOOST_CHECK(Utils::Vector3d::broadcast(1.) == box.length_inv());
    BOOST_CHECK(Utils::Vector3d::broadcast(0.5) == box.length_half());
  }

  /* setter */
  {
    auto box = BoxGeometry{};
    auto const box_l = Utils::Vector3d{1., 2., 3.};
    auto const box_l_inv =
        Utils::Vector3d{1. / box_l[0], 1. / box_l[1], 1. / box_l[2]};

    box.set_length(box_l);

    BOOST_CHECK(box.length() == box_l);
    BOOST_CHECK(box.length_inv() == box_l_inv);
    BOOST_CHECK(box.length_half() == 0.5 * box_l);
  }
}

namespace {
/**
 * @brief Build the deterministic set of per-axis separations swept by the
 * cuboid minimum-image tests.
 *
 * Covers, for a box length @p box_length: exact zero, small positive and
 * negative offsets, offsets just under +/- L/2, offsets in the wrap region
 * between L/2 and L, and offsets up to just under +/- 1.5 * L (the documented
 * validity domain of the branchless fold). The exact +/- L/2 ties are
 * deliberately excluded here; they are exercised by a separate tie test
 * because `rint` (ties-to-even) may legitimately pick the opposite image.
 */
std::vector<double> make_separations(double box_length) {
  std::vector<double> const fractions = {0.,    0.01, -0.01, 0.25, -0.25, 0.49,
                                         -0.49, 0.51, -0.51, 0.75, -0.75, 0.99,
                                         -0.99, 1.25, -1.25, 1.49, -1.49};
  std::vector<double> separations;
  separations.reserve(fractions.size());
  for (auto const fraction : fractions) {
    separations.push_back(fraction * box_length);
  }
  return separations;
}

/** @brief A deterministic pair of positions and the box periodicity mask. */
struct SweepPair {
  Utils::Vector3d a;
  Utils::Vector3d b;
};

/**
 * @brief Enumerate all deterministic position pairs of the sweep for a box.
 *
 * The initial point @p b is placed at a fixed offset (including negative
 * coordinates and coordinates outside the primary box), and the terminal
 * point @p a is obtained by adding every combination of per-axis separations.
 */
std::vector<SweepPair> make_sweep_pairs(Utils::Vector3d const &box_l) {
  auto const sep_x = make_separations(box_l[0u]);
  auto const sep_y = make_separations(box_l[1u]);
  auto const sep_z = make_separations(box_l[2u]);
  // base points: negative and outside the primary box
  std::array<Utils::Vector3d, 3u> const bases = {
      Utils::Vector3d{0.3 * box_l[0u], -1.7 * box_l[1u], 2.4 * box_l[2u]},
      Utils::Vector3d{-0.6 * box_l[0u], 0.2 * box_l[1u], -0.9 * box_l[2u]},
      Utils::Vector3d{1.1 * box_l[0u], 1.3 * box_l[1u], 0.05 * box_l[2u]}};
  std::vector<SweepPair> pairs;
  for (auto const &base : bases) {
    for (auto const dx : sep_x) {
      for (auto const dy : sep_y) {
        for (auto const dz : sep_z) {
          auto const b = base;
          auto const a = Utils::Vector3d{b[0u] + dx, b[1u] + dy, b[2u] + dz};
          pairs.push_back(SweepPair{a, b});
        }
      }
    }
  }
  return pairs;
}
} // namespace

BOOST_AUTO_TEST_CASE(cuboid_minimum_image_single_pair_test) {
  auto const box_l = Utils::Vector3d{2.5, 4.0, 7.25};
  auto const pairs = make_sweep_pairs(box_l);

  for (unsigned mask = 0u; mask < 8u; ++mask) {
    auto box = BoxGeometry{};
    box.set_type(BoxType::CUBOID);
    box.set_length(box_l);
    box.set_periodic(0u, (mask & 1u) != 0u);
    box.set_periodic(1u, (mask & 2u) != 0u);
    box.set_periodic(2u, (mask & 4u) != 0u);

    auto const fold = box.cuboid_minimum_image();

    for (auto const &pair : pairs) {
      auto const &a = pair.a;
      auto const &b = pair.b;
      auto const reference = box.get_mi_vector(a, b);

      // dist2 must match get_mi_vector(...).norm2() bitwise.
      BOOST_CHECK_EQUAL(fold.dist2(a, b), reference.norm2());
    }
  }
}

BOOST_AUTO_TEST_CASE(cuboid_minimum_image_half_box_tie_test) {
  auto const box_l = Utils::Vector3d{2.5, 4.0, 7.25};

  // periodic in all directions: separation exactly L/2 on each axis.
  auto box = BoxGeometry{};
  box.set_type(BoxType::CUBOID);
  box.set_length(box_l);

  auto const fold = box.cuboid_minimum_image();

  // Place b at the origin and a at exactly L/2 so the per-axis separation is
  // bitwise-exactly half a box length.
  auto const b = Utils::Vector3d{0., 0., 0.};
  auto const a = box.length_half();

  auto const reference = box.get_mi_vector(a, b);
  // At the exact tie, only the distance is well-defined (both images are
  // equidistant); rint may pick the opposite sign of the standard path.
  BOOST_CHECK_EQUAL(fold.dist2(a, b), reference.norm2());
  // Each periodic axis must contribute exactly (L/2)^2.
  auto const expected = box.length_half().norm2();
  BOOST_CHECK_EQUAL(fold.dist2(a, b), expected);
}

BOOST_AUTO_TEST_CASE(cuboid_minimum_image_batch_test) {
  auto const box_l = Utils::Vector3d{2.5, 4.0, 7.25};
  auto const pairs = make_sweep_pairs(box_l);

  // all periodic: exercises folding on every axis.
  auto box = BoxGeometry{};
  box.set_type(BoxType::CUBOID);
  box.set_length(box_l);
  auto const fold = box.cuboid_minimum_image();

  // A single reference point; the neighbour SoA arrays hold the terminal
  // points of the sweep. Chunk into tiles, including a non-full final tile
  // and a tile of size one.
  Utils::Vector3d const reference_point =
      Utils::Vector3d{0.13 * box_l[0u], -0.87 * box_l[1u], 1.9 * box_l[2u]};

  std::vector<double> sx, sy, sz;
  sx.reserve(pairs.size());
  sy.reserve(pairs.size());
  sz.reserve(pairs.size());
  for (auto const &pair : pairs) {
    sx.push_back(pair.a[0u]);
    sy.push_back(pair.a[1u]);
    sz.push_back(pair.a[2u]);
  }

  auto const n = static_cast<int>(pairs.size());
  std::vector<int> const tile_sizes = {1, 7, 8, 16, 15};

  int offset = 0;
  std::size_t tile_index = 0u;
  while (offset < n) {
    auto const capacity = tile_sizes[tile_index % tile_sizes.size()];
    ++tile_index;
    auto const m = std::min(capacity, n - offset); // last tile is partial
    std::vector<double> dx0(static_cast<std::size_t>(capacity)),
        dx1(static_cast<std::size_t>(capacity)),
        dx2(static_cast<std::size_t>(capacity)),
        dsq(static_cast<std::size_t>(capacity));

    fold.batch_vector_dist2(reference_point[0u], reference_point[1u],
                            reference_point[2u], m, sx.data() + offset,
                            sy.data() + offset, sz.data() + offset, dx0.data(),
                            dx1.data(), dx2.data(), dsq.data());

    for (int t = 0; t < m; ++t) {
      auto const idx =
          static_cast<std::size_t>(offset) + static_cast<std::size_t>(t);
      auto const neighbour = Utils::Vector3d{sx[idx], sy[idx], sz[idx]};
      // get_mi_vector is the scalar oracle: for cuboid boxes it runs the
      // same per-axis get_mi_coord_masked fold as the batch path.
      auto const single = box.get_mi_vector(reference_point, neighbour);
      auto const single_dist2 = fold.dist2(reference_point, neighbour);
      BOOST_CHECK_EQUAL(dx0[static_cast<std::size_t>(t)], single[0u]);
      BOOST_CHECK_EQUAL(dx1[static_cast<std::size_t>(t)], single[1u]);
      BOOST_CHECK_EQUAL(dx2[static_cast<std::size_t>(t)], single[2u]);
      BOOST_CHECK_EQUAL(dsq[static_cast<std::size_t>(t)], single_dist2);
    }

    offset += m;
  }
}

BOOST_AUTO_TEST_CASE(cuboid_minimum_image_non_periodic_axis_test) {
  auto const box_l = Utils::Vector3d{2.5, 4.0, 7.25};

  // periodic in x and z, non-periodic in y.
  auto box = BoxGeometry{};
  box.set_type(BoxType::CUBOID);
  box.set_length(box_l);
  box.set_periodic(0u, true);
  box.set_periodic(1u, false);
  box.set_periodic(2u, true);
  auto const fold = box.cuboid_minimum_image();

  // Separations larger than L on every axis.
  auto const b = Utils::Vector3d{0.1, 0.2, 0.3};
  auto const raw =
      Utils::Vector3d{1.3 * box_l[0u], 1.3 * box_l[1u], 1.3 * box_l[2u]};
  auto const a = b + raw;

  auto const reference = box.get_mi_vector(a, b);

  // The non-periodic y axis must not fold: component equals the raw
  // separation.
  BOOST_CHECK_EQUAL(reference[1u], raw[1u]);
  // The periodic axes must fold (magnitude reduced below L/2).
  BOOST_CHECK_LT(std::abs(reference[0u]), 0.5 * box_l[0u]);
  BOOST_CHECK_LT(std::abs(reference[2u]), 0.5 * box_l[2u]);

  // And the hoisted fold must agree with the standard box path bitwise.
  BOOST_CHECK_EQUAL(fold.dist2(a, b), reference.norm2());
}
