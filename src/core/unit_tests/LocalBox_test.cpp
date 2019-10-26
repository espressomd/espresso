/*
 * Copyright (C) 2019 The ESPResSo project
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
#include <boost/range/algorithm/equal.hpp>
#include <boost/test/unit_test.hpp>

#include <LocalBox.hpp>

#include <limits>

/**
 * @brief Check that the box corners and side length
 *        agree.
 */
template <class T> void check_length(LocalBox<T> box) {
  auto const expected = box.my_right() - box.my_left();
  auto const result = box.length();

  BOOST_CHECK_SMALL((result - expected).norm2(),
                    std::numeric_limits<T>::epsilon());
}

BOOST_AUTO_TEST_CASE(constructors) {
  /* default */
  {
    auto const box = LocalBox<float>();
    BOOST_CHECK_EQUAL(box.my_left().norm2(), 0.f);
    check_length(box);
  }

  /* corner + length */
  {
    Utils::Vector<double, 3> const lower_corner = {1., 2., 3.};
    Utils::Vector<double, 3> const local_box_length = {4., 5., 6.};
    Utils::Array<int, 6> const boundaries = {-1, 0, 1, 1, 0, -1};

    auto const box =
        LocalBox<double>(lower_corner, local_box_length, boundaries);

    BOOST_CHECK(box.my_left() == lower_corner);
    BOOST_CHECK(box.length() == local_box_length);
    BOOST_CHECK(boost::equal(boundaries, box.boundary()));
    check_length(box);
  }
}