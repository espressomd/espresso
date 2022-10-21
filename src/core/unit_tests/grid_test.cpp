/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

#include <utils/Vector.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

template <class T> static auto epsilon = std::numeric_limits<T>::epsilon();

BOOST_AUTO_TEST_CASE(get_mi_coord_test) {
  using detail::get_mi_coord;

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
  using detail::get_mi_coord;
  using Utils::Vector3d;

  Vector3d box_l = {1., 2., 3.};
  BoxGeometry box;
  box.set_length(box_l);
  box.set_periodic(1, false);

  auto const a = Vector3d{1.1, 12.2, -13.4};
  auto const b = Vector3d{-0.9, 8.8, 21.1};

  auto const result = box.get_mi_vector(a, b);

  for (int i = 0; i < 3; i++) {
    auto const expected = get_mi_coord(a[i], b[i], box_l[i], box.periodic(i));

    BOOST_CHECK_SMALL(std::abs(expected - result[i]), epsilon<double>);
  }
}

BOOST_AUTO_TEST_CASE(lees_edwards_mi_vector) {
  using detail::get_mi_coord;
  using Utils::Vector3d;

  Vector3d box_l = {5., 2., 3.};
  BoxGeometry box;
  box.set_length(box_l);
  box.set_periodic(1, false);
  box.set_type(BoxType::LEES_EDWARDS);
  LeesEdwardsBC le{0., 0., 2, 0};
  box.set_lees_edwards_bc(le);

  // No pos offset -> behave like normal get_mi_vector
  {
    auto const a = Vector3d{5.1, 12.2, -13.4};
    auto const b = Vector3d{-0.9, 8.8, 21.1};

    auto const result = box.get_mi_vector(a, b);

    for (int i = 0; i < 3; i++) {
      auto const expected = get_mi_coord(a[i], b[i], box_l[i], box.periodic(i));
      BOOST_CHECK_SMALL(std::abs(expected - result[i]), 5. * epsilon<double>);
    }
  }

  // LE pos offset >0 but distance < half box length. Behave like normal
  // get_mi_vector
  le.pos_offset = 1.;
  box.set_lees_edwards_bc(le);
  {
    auto const a = Vector3d{1.1, 12.2, -13.4};
    auto const b = Vector3d{-0.9, 8.8, 21.1};

    auto const result = box.get_mi_vector(a, b);

    for (int i = 0; i < 3; i++) {
      auto expected = get_mi_coord(a[i], b[i], box_l[i], box.periodic(i));
      if (i == le.shear_direction) {
        expected -= le.pos_offset = 1.;
      }
      BOOST_CHECK_SMALL(std::abs(expected - result[i]), 5. * epsilon<double>);
    }
  }
  // LE pos offset and distance > box in shear plane normal direction
  {
    auto const a = Vector3d{1.1, 12.2, -13.};
    auto const b = Vector3d{-11., 8.8, -13.};
    auto const le_jumps = std::round((a[0] - b[0]) / box.length()[0]);

    auto const result = box.get_mi_vector(a, b);

    auto const expected = Vector3d{{
        std::fmod(a[0] - b[0], box.length()[0]), a[1] - b[1],
        a[2] - b[2] + le_jumps * le.pos_offset -
            box.length()[2] // Manually apply minimum image convention
    }};
    for (int i : {0, 1, 2}) {
      BOOST_CHECK_CLOSE(expected[i], result[i], 100. * epsilon<double>);
    }
  }

  // Test a case where coordinate different + LE offset in shift dir > box_l/2
  box.set_type(BoxType::LEES_EDWARDS);
  box.set_periodic(0, true);
  box.set_periodic(1, true);
  box.set_periodic(2, true);
  box.set_length(Vector3d{5., 5., 5.});
  box.set_lees_edwards_bc(LeesEdwardsBC{2.98, 0., 0, 1});
  auto const result =
      box.get_mi_vector(Vector3d{2.5, 1., 2.5}, Vector3d{4.52, 4., 2.5});
  BOOST_CHECK_SMALL(std::fabs(result[0]), box.length_half()[0]);
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
  Utils::Vector3d pos{5., 6., 7.};
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

BOOST_AUTO_TEST_CASE(fold_position_test) {
  auto const box_l = Utils::Vector3d{2., 4., 6.};
  auto box = BoxGeometry();
  box.set_length(box_l);
  box.set_periodic(0, true);
  box.set_periodic(1, true);
  box.set_periodic(2, false);

  /* Wrapped */
  {
    Utils::Vector3d pos{-1.9, 4.1, 7.};
    Utils::Vector3i img{0, 1, 1};
    Utils::Vector3d const expected_pos{0.1, 0.1, 7.};
    Utils::Vector3i const expected_img{-1, 2, 1};

    fold_position(pos, img, box);

    BOOST_CHECK_SMALL((pos - expected_pos).norm(), 3 * epsilon<double>);
    BOOST_CHECK_EQUAL((img - expected_img).norm2(), 0);
  }

  /* Not wrapped */
  {
    Utils::Vector3d pos{1., 3., 7.};
    Utils::Vector3i img{0, -1, -1};
    Utils::Vector3d const expected_pos{1., 3., 7.};
    Utils::Vector3i const expected_img{0, -1, -1};

    fold_position(pos, img, box);

    BOOST_CHECK_SMALL((pos - expected_pos).norm(), 3 * epsilon<double>);
    BOOST_CHECK_EQUAL((img - expected_img).norm2(), 0);
  }
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
