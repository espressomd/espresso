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

#define BOOST_TEST_MODULE field index range
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "../src/utils/field_index_range.hpp"

#include <core/DataTypes.h>

#include <cstddef>
#include <stdexcept>
#include <string>

BOOST_AUTO_TEST_CASE(index_range_of_the_padded_field) {
  using walberla::uint_t;
  using walberla::detail::check_field_index_range;
  auto constexpr limit = uint_t{1} << 31u;
  auto const check_range = [](std::size_t n_components, uint_t cells_x,
                              uint_t cells_y, uint_t cells_z,
                              uint_t n_ghost_layers) {
    check_field_index_range("test field", "test kernels", n_components, cells_x,
                            cells_y, cells_z, n_ghost_layers);
  };
  // one component, no ghost layers: the field has pad(x) * y * z elements,
  // where the x extent is padded to a multiple of 16 elements (64 bytes of
  // floats); the padding makes the difference in the first pair
  BOOST_CHECK_NO_THROW(check_range(1u, limit - 16u, 1u, 1u, 0u));
  BOOST_CHECK_THROW(check_range(1u, limit - 15u, 1u, 1u, 0u),
                    std::runtime_error);
  // only the x extent is padded
  BOOST_CHECK_NO_THROW(check_range(1u, 16u, (limit / 16u) - 1u, 1u, 0u));
  BOOST_CHECK_THROW(check_range(1u, 16u, limit / 16u, 1u, 0u),
                    std::runtime_error);
  // the ghost layers are part of the size, on both sides of each extent:
  // extents 16 * 2^14 * (2^13 - 1) < 2^31 <= 16 * 2^14 * 2^13
  BOOST_CHECK_NO_THROW(
      check_range(1u, 14u, (1u << 14u) - 2u, (1u << 13u) - 3u, 1u));
  BOOST_CHECK_THROW(
      check_range(1u, 14u, (1u << 14u) - 2u, (1u << 13u) - 2u, 1u),
      std::runtime_error);
  // two ghost layers: x extent 8 + 4 = 12 is padded to 16
  BOOST_CHECK_NO_THROW(
      check_range(1u, 8u, (1u << 14u) - 4u, (1u << 13u) - 5u, 2u));
  BOOST_CHECK_THROW(check_range(1u, 8u, (1u << 14u) - 4u, (1u << 13u) - 4u, 2u),
                    std::runtime_error);
  // the components are part of the size: 19 * 16 * 4096 * 1724 < 2^31
  BOOST_CHECK_NO_THROW(check_range(19u, 14u, 4094u, 1722u, 1u));
  BOOST_CHECK_THROW(check_range(19u, 14u, 4094u, 1723u, 1u),
                    std::runtime_error);
  // cubic blocks of 19 components with one ghost layer: for 479^3 the x extent
  // 481 is padded to 496, which pushes the field over the limit
  BOOST_CHECK_NO_THROW(check_range(19u, 478u, 478u, 478u, 1u));
  BOOST_CHECK_THROW(check_range(19u, 479u, 479u, 479u, 1u), std::runtime_error);
  // non-cubic blocks of 19 components with one ghost layer: x extent 64, y
  // extent 702
  BOOST_CHECK_NO_THROW(check_range(19u, 62u, 700u, 2513u, 1u));
  BOOST_CHECK_THROW(check_range(19u, 62u, 700u, 2514u, 1u), std::runtime_error);
  // the size computation doesn't overflow
  BOOST_CHECK_THROW(check_range(19u, uint_t{1} << 62u, uint_t{1} << 62u,
                                uint_t{1} << 62u, 1u),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(index_range_of_the_flux_field) {
  using walberla::uint_t;
  using walberla::detail::check_field_index_range;
  auto const check_flux = [](uint_t cells_x, uint_t cells_y, uint_t cells_z,
                             uint_t n_ghost_layers) {
    check_field_index_range("flux field", "EK kernels", 13u, cells_x, cells_y,
                            cells_z, n_ghost_layers);
  };
  // cubic blocks with one ghost layer: for 542^3 the x extent 544 is padded
  // to 544 and 13 * 544 * 544^2 < 2^31, for 543^3 the x extent 545 is padded
  // to 560 and 13 * 560 * 545^2 > 2^31
  BOOST_CHECK_NO_THROW(check_flux(542u, 542u, 542u, 1u));
  BOOST_CHECK_THROW(check_flux(543u, 543u, 543u, 1u), std::runtime_error);
  // two ghost layers: the extents are 540 + 4 = 544 and 541 + 4 = 545
  BOOST_CHECK_NO_THROW(check_flux(540u, 540u, 540u, 2u));
  BOOST_CHECK_THROW(check_flux(541u, 541u, 541u, 2u), std::runtime_error);
  // the same blocks with a field of one component (the density field)
  BOOST_CHECK_NO_THROW(check_field_index_range("density field", "EK kernels",
                                               1u, 543u, 543u, 543u, 1u));
}

BOOST_AUTO_TEST_CASE(error_message) {
  using walberla::detail::check_field_index_range;
  auto const check_message = [](std::runtime_error const &error) {
    std::string const message(error.what());
    return message.find("The flux field of a block (13 components in "
                        "600x601x602 cells with 2 ghost layers)") == 0u and
           message.find("2^31") != std::string::npos and
           message.find("of the EK kernels") != std::string::npos;
  };
  BOOST_CHECK_EXCEPTION(check_field_index_range("flux field", "EK kernels", 13u,
                                                600u, 601u, 602u, 2u),
                        std::runtime_error, check_message);
}
