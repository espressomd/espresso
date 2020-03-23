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
#include <boost/test/unit_test.hpp>

#include "grid.hpp"

#include <cmath>
#include <limits>

template <class T> auto const epsilon = std::numeric_limits<T>::epsilon();

BOOST_AUTO_TEST_CASE(get_mi_coord_test) {
  auto const box_l = 3.1415;

  /* Non-periodic */
  {{auto const a = 1.;
  auto const b = 2.;

  BOOST_CHECK_EQUAL(get_mi_coord(a, b, box_l, /* periodic */ false), a - b);
}

{
  auto const a = 1.;
  auto const b = 3.;

  BOOST_CHECK_EQUAL(get_mi_coord(a, b, box_l, /* periodic */ false), a - b);
}
}

/* Regular distance */
{
  auto const a = -0.5;
  auto const b = +1.0;

  BOOST_CHECK_EQUAL(get_mi_coord(a, b, box_l, /* periodic */ true), a - b);
  BOOST_CHECK_EQUAL(get_mi_coord(b, a, box_l, /* periodic */ true), b - a);
}

/* Wrapped */
{
  auto const a = 1.;
  auto const b = 3.;

  BOOST_CHECK_SMALL(std::abs(get_mi_coord(a, b, box_l, /* periodic */ true) -
                             (a - b) - box_l),
                    epsilon<double>);
  BOOST_CHECK_SMALL(std::abs(get_mi_coord(b, a, box_l, /* periodic */ true) -
                             (b - a) + box_l),
                    epsilon<double>);
}

/* Corner cases */
{
  {
    auto const a = 0.4;
    auto const b = a + 0.5 * box_l;

    BOOST_CHECK_SMALL(std::abs(get_mi_coord(a, b, box_l, /* periodic */ true) -
                               (a - b) - box_l),
                      epsilon<double>);
    BOOST_CHECK_SMALL(std::abs(get_mi_coord(b, a, box_l, /* periodic */ true) -
                               (b - a) + box_l),
                      epsilon<double>);
  }

  {
    auto const a = 0.4;
    auto const b = std::nextafter(a + 0.5 * box_l, box_l);

    BOOST_CHECK_SMALL(std::abs(get_mi_coord(a, b, box_l, /* periodic */ true) -
                               (a - b) - box_l),
                      epsilon<double>);
    BOOST_CHECK_SMALL(std::abs(get_mi_coord(b, a, box_l, /* periodic */ true) -
                               (b - a) + box_l),
                      epsilon<double>);
  }

  {
    auto const a = 0.4;
    auto const b = std::nextafter(a + 0.5 * box_l, 0.);

    BOOST_CHECK_SMALL(
        std::abs(get_mi_coord(a, b, box_l, /* periodic */ true) - (a - b)),
        epsilon<double>);
    BOOST_CHECK_SMALL(
        std::abs(get_mi_coord(b, a, box_l, /* periodic */ true) - (b - a)),
        epsilon<double>);
  }
}
}

BOOST_AUTO_TEST_CASE(get_mi_vector_test) {
  using Utils::Vector3d;

  Vector3d box_l = {1., 2., 3.};
  BoxGeometry box;
  box.set_length(box_l);
  box.set_periodic(1, false);

  auto const a = Vector3d{1.1, 12.2, -13.4};
  auto const b = Vector3d{-0.9, 8.8, 21.1};

  auto const result = get_mi_vector(a, b, box);

  for (int i = 0; i < 3; i++) {
    auto const expected = get_mi_coord(a[i], b[i], box_l[i], box.periodic(i));

    BOOST_CHECK_SMALL(std::abs(expected - result[i]), epsilon<double>);
  }
}

BOOST_AUTO_TEST_CASE(image_shift_test) {
  Utils::Vector3i img{1, -2, 3};
  Utils::Vector3d box{1., 2., 3.};

  auto const result = image_shift(img, box);
  auto const expected =
      Utils::Vector3d{img[0] * box[0], img[1] * box[1], img[2] * box[2]};

  BOOST_CHECK_SMALL((result - expected).norm(), epsilon<double>);
}

BOOST_AUTO_TEST_CASE(unfolded_position_test) {
  Utils::Vector3d pos{5., 6, 7.};
  Utils::Vector3i img{1, -2, 3};
  Utils::Vector3d box{1., 2., 3.};

  auto expected = pos + image_shift(img, box);
  auto result = unfolded_position(pos, img, box);

  BOOST_CHECK_SMALL((result - expected).norm(), epsilon<double>);
}

BOOST_AUTO_TEST_CASE(fold_coordinate_test) {
  BOOST_CHECK_THROW(fold_coordinate(0., std::numeric_limits<int>::max(), 1.),
                    std::runtime_error);
  BOOST_CHECK_THROW(fold_coordinate(0., std::numeric_limits<int>::min(), 1.),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(regular_decomposition_test) {
  auto const eps = std::numeric_limits<double>::epsilon();

  auto const box_l = Utils::Vector3d{10, 20, 30};
  auto box = BoxGeometry();
  box.set_length(box_l);
  auto const node_grid = Utils::Vector3i{1, 2, 3};

  /* check length */
  {
    auto const result = regular_decomposition(box, {0, 0, 0}, node_grid);
    auto const local_box_l = result.length();

    BOOST_CHECK_CLOSE(box_l[0], local_box_l[0] * node_grid[0], 100. * eps);
    BOOST_CHECK_CLOSE(box_l[1], local_box_l[1] * node_grid[1], 100. * eps);
    BOOST_CHECK_CLOSE(box_l[2], local_box_l[2] * node_grid[2], 100. * eps);
  }

  /* check corners */
  {
    Utils::Vector3i node_pos;
    for (node_pos[0] = 0; node_pos[0] < node_grid[0]; node_pos[0]++)
      for (node_pos[1] = 0; node_pos[1] < node_grid[1]; node_pos[1]++)
        for (node_pos[2] = 0; node_pos[2] < node_grid[2]; node_pos[2]++) {
          auto const result = regular_decomposition(box, node_pos, node_grid);
          auto const local_box_l = result.length();
          auto const lower_corner = result.my_left();

          BOOST_CHECK_CLOSE(lower_corner[0], local_box_l[0] * node_pos[0],
                            100. * eps);
          BOOST_CHECK_CLOSE(lower_corner[1], local_box_l[1] * node_pos[1],
                            100. * eps);
          BOOST_CHECK_CLOSE(lower_corner[2], local_box_l[2] * node_pos[2],
                            100. * eps);
        }
  }
}
