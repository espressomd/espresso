/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE coordinate_transformation test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/coordinate_transformation.hpp>
#include <utils/math/vec_rotate.hpp>

#include <cmath>

using Utils::Vector3d;

BOOST_AUTO_TEST_CASE(cartesian_to_cylinder_test) {
  Vector3d const cart_coord{{1.0, 3.3, 2.0}};
  auto const transformed_x = transform_coordinate_cartesian_to_cylinder(
      cart_coord, Vector3d{{1, 0, 0}});
  auto const transformed_y = transform_coordinate_cartesian_to_cylinder(
      cart_coord, Vector3d{{0, 1, 0}});
  auto const transformed_z = transform_coordinate_cartesian_to_cylinder(
      cart_coord, Vector3d{{0, 0, 1}});
  // For x as the symmetry axis we rotate the cartesian coordinates around the
  // y-axis by -pi/2.
  auto const expected_x = transform_coordinate_cartesian_to_cylinder(
      vec_rotate(Vector3d{{0.0, 1.0, 0.0}}, -Utils::pi() / 2.0, cart_coord),
      Vector3d{{0, 0, 1}});
  // For y as the symmetry axis we rotate the cartesian coordinates around the
  // x-axis by pi/2.
  auto const expected_y = transform_coordinate_cartesian_to_cylinder(
      vec_rotate(Vector3d{{1.0, 0.0, 0.0}}, Utils::pi() / 2.0, cart_coord),
      Vector3d{{0, 0, 1}});
  auto const expected_z = Vector3d{
      {std::sqrt(cart_coord[0] * cart_coord[0] + cart_coord[1] * cart_coord[1]),
       std::atan2(cart_coord[1], cart_coord[0]), cart_coord[2]}};

  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK(transformed_x[i] == expected_x[i]);
    BOOST_CHECK(transformed_y[i] == expected_y[i]);
    BOOST_CHECK(transformed_z[i] == expected_z[i]);
  }
}

BOOST_AUTO_TEST_CASE(cylinder_to_cartesian_test) {
  Vector3d const cylinder_coord{{1.2, 3.123, 42.0}};
  auto const transformed_x = transform_coordinate_cylinder_to_cartesian(
      cylinder_coord, Vector3d{{1, 0, 0}});
  auto const transformed_y = transform_coordinate_cylinder_to_cartesian(
      cylinder_coord, Vector3d{{0, 1, 0}});
  auto const transformed_z = transform_coordinate_cylinder_to_cartesian(
      cylinder_coord, Vector3d{{0, 0, 1}});
  // We transform from cylinder zu cartesian and have to rotate back. See test
  // cartesian_to_cylinder_test.
  auto const expected_x =
      vec_rotate(Vector3d{{0.0, 1.0, 0.0}}, Utils::pi() / 2.0,
                 transform_coordinate_cylinder_to_cartesian(
                     cylinder_coord, Vector3d{{0, 0, 1}}));
  auto const expected_y =
      vec_rotate(Vector3d{{1.0, 0.0, 0.0}}, -Utils::pi() / 2.0,
                 transform_coordinate_cylinder_to_cartesian(
                     cylinder_coord, Vector3d{{0, 0, 1}}));
  // x = r * cos(phi); y = r * sin(phi); z = z
  auto const expected_z = Vector3d{
      {cylinder_coord[0] * std::cos(cylinder_coord[1]),
       cylinder_coord[0] * std::sin(cylinder_coord[1]), cylinder_coord[2]}};
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK(transformed_x[i] == expected_x[i]);
    BOOST_CHECK(transformed_y[i] == expected_y[i]);
    BOOST_CHECK(transformed_z[i] == expected_z[i]);
  }
}

BOOST_AUTO_TEST_CASE(vector_cart_to_cyl_test) {
  constexpr auto eps = 1e-13;
  Vector3d const pos{{1.1, 2.2, 3.3}};
  auto const axis = (Vector3d{{4.4, 5.5, 6.6}}).normalized();
  Vector3d const vec{{7.7, 8.8, 9.9}};

  auto const vec_cyl = transform_vector_cartesian_to_cylinder(vec, axis, pos);

  // cylindrical basis vectors at pos
  auto const e_z = axis;
  auto const e_r = (pos - (pos * axis) * axis).normalized();
  auto const e_phi = Utils::vector_product(e_z, e_r);

  BOOST_CHECK_SMALL(vec_cyl[0] - vec * e_r, eps);
  BOOST_CHECK_SMALL(vec_cyl[1] - vec * e_phi, eps);
  BOOST_CHECK_SMALL(vec_cyl[2] - vec * e_z, eps);
}
