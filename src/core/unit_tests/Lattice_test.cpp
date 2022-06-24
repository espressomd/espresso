/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE Lattice class tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "grid_based_algorithms/lattice.hpp"

#include <utils/Vector.hpp>

#include <cstddef>
#include <limits>
#include <stdexcept>

BOOST_AUTO_TEST_CASE(test_basic_lattice) {
  // create a lattice for the second domain of a 2x1x1 partition of the box
  auto const halo_size = Lattice::index_t{1};
  auto const agrid = 0.5;
  auto const offset = 0.5;
  auto const box_length = Utils::Vector3d{{6., 6., 6.}};
  auto const local_box = Utils::Vector3d{{3., 6., 6.}};
  auto const node_pos = Utils::Vector3i{{1, 0, 0}};
  auto const node_grid = Utils::Vector3i{{2, 1, 1}};
  Lattice lattice(agrid, offset, halo_size, local_box, box_length, box_length,
                  node_pos, node_grid);

  // check struct members
  BOOST_CHECK_EQUAL(lattice.halo_size, halo_size);
  BOOST_CHECK_EQUAL(lattice.agrid, agrid);
  BOOST_CHECK_EQUAL(lattice.offset, offset);
  BOOST_CHECK_EQUAL(lattice.halo_grid_volume, (6 + 2) * (12 + 2) * (12 + 2));
  auto const elementwise = boost::test_tools::per_element();
  auto const ref_grid = Utils::Vector3i{{6, 12, 12}};
  auto const ref_global_grid = Utils::hadamard_product(node_grid, ref_grid);
  auto const local_index_offset = Utils::hadamard_product(node_pos, ref_grid);
  BOOST_TEST(lattice.local_box == local_box, elementwise);
  BOOST_TEST(lattice.node_grid == node_grid, elementwise);
  BOOST_TEST(lattice.grid == ref_grid, elementwise);
  BOOST_TEST(lattice.global_grid == ref_global_grid, elementwise);
  BOOST_TEST(lattice.local_index_offset == local_index_offset, elementwise);

  // check methods
  BOOST_CHECK(lattice.is_local({11, 11, 11}));
  BOOST_CHECK(lattice.is_local({6, 11, 11}));
  BOOST_CHECK(!lattice.is_local({5, 11, 11}));
  BOOST_CHECK(!lattice.is_local({12, 12, 12}));
  BOOST_CHECK(!lattice.is_local({0, 0, 0}));
  auto const global_index = Utils::Vector3i{{11, 11, 11}};
  auto const local_index = Utils::Vector3i{{6, 12, 12}};
  BOOST_TEST(lattice.local_index(global_index) == local_index, elementwise);
}

BOOST_AUTO_TEST_CASE(test_map_position_to_lattice) {
  using boost::test_tools::per_element;
  auto const halo_size = Lattice::index_t{1};
  auto const agrid = 1.0;
  auto const offset = 0.5;
  auto const box_l = Utils::Vector3d{{6., 6., 6.}};
  auto const local_box = Utils::Vector3d{{6., 6., 6.}};
  auto const node_pos = Utils::Vector3i{{0, 0, 0}};
  auto const node_grid = Utils::Vector3i{{1, 1, 1}};
  Lattice lattice(agrid, offset, halo_size, local_box, box_l, box_l, node_pos,
                  node_grid);

  // check methods
  auto const slice_x = 6u + 2u;
  auto const slice_xy = slice_x * slice_x;
  auto const slice_xyz = 2u * 6u * 6u;
  Utils::Vector<std::size_t, 8> const origin_index = {
      0u,        1u,
      slice_x,   slice_x + 1u,
      slice_xy,  slice_xy + 1u,
      slice_xyz, slice_xyz + 1u};
  auto const delta1_ref = Utils::Vector6d{{.5, .5, .5, .5, .5, .5}};
  auto const delta2_ref = Utils::Vector6d{{1., 1., 1., 0., 0., 0.}};
  Utils::Vector<std::size_t, 8> node_index1;
  Utils::Vector<std::size_t, 8> node_index2;
  Utils::Vector<std::size_t, 8> idx;
  Utils::Vector6d delta1;
  Utils::Vector6d delta2;
  Utils::Vector6d dx;

  // check inside local domain (edge cases)
  auto const my_origin = Utils::Vector3d::broadcast(0.);
  auto const my_lb_left = Utils::Vector3d::broadcast(-offset);
  auto const my_lb_right = Utils::Vector3d::broadcast(offset - 1e-12) + box_l;
  lattice.map_position_to_lattice(my_origin, node_index1, delta1);
  lattice.map_position_to_lattice(my_lb_left, node_index2, delta2);
  lattice.map_position_to_lattice(my_lb_right, idx, dx);
  BOOST_TEST(node_index1 == origin_index, per_element());
  BOOST_TEST(node_index2 == origin_index, per_element());
  BOOST_TEST(delta1 == delta1_ref, per_element());
  BOOST_TEST(delta2 == delta2_ref, per_element());

  // check almost inside local domain
  auto constexpr epsilon = std::numeric_limits<double>::epsilon();
  if (epsilon != epsilon / 2.) { // check for machine precision
    auto const outside = Utils::Vector3d::broadcast(-offset - epsilon / 2.);
    lattice.map_position_to_lattice(outside, node_index2, delta2);
    BOOST_TEST(node_index2 == origin_index, per_element());
  }

  // check outside local domain
  BOOST_CHECK_THROW(lattice.map_position_to_lattice({-2., -2., -2.}, idx, dx),
                    std::runtime_error);
  BOOST_CHECK_THROW(lattice.map_position_to_lattice({6.5, 6.5, 6.5}, idx, dx),
                    std::runtime_error);
}
