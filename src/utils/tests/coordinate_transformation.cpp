/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE coordinate_transformation test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/coordinate_transformation.hpp>
#include <utils/math/vec_rotate.hpp>

#include <cmath>
#include <random>

using Utils::Vector3d;

BOOST_AUTO_TEST_CASE(cartesian_to_cylinder_test) {
  constexpr auto eps = 1e-14;
  auto const pos = Vector3d{{1.0, 3.3, 2.0}};
  auto const cyl = transform_coordinate_cartesian_to_cylinder(pos);
  BOOST_CHECK_SMALL(cyl[0] - std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]), eps);
  BOOST_CHECK_SMALL(cyl[1] - std::atan2(pos[1], pos[0]), eps);
  BOOST_CHECK_SMALL(cyl[2] - pos[2], eps);
}

BOOST_AUTO_TEST_CASE(basis_transform_test) {
  constexpr auto eps = 1e-14;
  Vector3d const b_x{{1, 0, 0}};
  Vector3d const b_y{{0, 1, 0}};
  Vector3d const b_z{{0, 0, 1}};
  // identity transform
  Vector3d const v{{1, 2, 3}};
  Vector3d const v_identity_transform = Utils::basis_change(b_x, b_y, b_z, v);
  // identity transform (swap both the vector and coordinate system)
  Vector3d const v_swap_coord_transform =
      Utils::basis_change(b_z, b_y, b_x, {{v[2], v[1], v[0]}});
  // non-trivial transform
  Vector3d const v1 = Vector3d{{2, 2, 2}}.normalized();
  Vector3d const v2 = Vector3d{{3, 3, -6}}.normalized();
  Vector3d const v3 = Utils::vector_product(v1, v2).normalized();
  Vector3d const v4 = basis_change(v1, v2, v3, 0.1 * v1 + 0.2 * v2 - 0.3 * v3);
  Vector3d const v4_expected = Vector3d{{0.1, 0.2, -0.3}};
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_SMALL(v_identity_transform[i] - v[i], eps);
    BOOST_CHECK_SMALL(v_swap_coord_transform[i] - v[i], eps);
    BOOST_CHECK_SMALL(v4[i] - v4_expected[i], eps);
  }
}

BOOST_AUTO_TEST_CASE(
    transform_coordinate_cartesian_to_cylinder_base_change_test) {
  constexpr auto eps = 1e-14;
  Vector3d const v1{{1, 3, 4}};
  Vector3d const v2{{-3.0, 7, 2}};
  Vector3d const axis = Utils::vector_product(v1, v2).normalized();
  Vector3d const v3 =
      Utils::transform_coordinate_cartesian_to_cylinder(v1, axis, v1);
  Vector3d const v3_ref{{v1.norm(), 0, 0}};
  auto const angle_v1_v2 = Utils::angle_between(v1, v2);
  auto const v4 = Utils::transform_coordinate_cartesian_to_cylinder(
      v2 + 2 * axis, axis, v1);
  Vector3d v4_ref{{v2.norm(), angle_v1_v2, 2}};
  auto const v5 = Utils::transform_coordinate_cartesian_to_cylinder(
      v1 + 2 * axis, axis, v2);
  Vector3d v5_ref{{v1.norm(), -angle_v1_v2, 2}};
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_SMALL(v3[i] - v3_ref[i], eps);
    BOOST_CHECK_SMALL(v4[i] - v4_ref[i], eps);
    BOOST_CHECK_SMALL(v5[i] - v5_ref[i], eps);
  }
}

BOOST_AUTO_TEST_CASE(cartesian_to_cylinder_with_axis_and_orientation_test) {
  constexpr auto eps = 1e-14;
  // tilted orthogonal basis
  auto const y = (Vector3d{{0, 1, -1}}).normalize();
  auto const z = (Vector3d{{1, 1, 1}}).normalize();
  auto const x = Utils::vector_product(y, z);

  // check transformation with orientation (phi is random for r=0)
  {
    auto const x_cyl = transform_coordinate_cartesian_to_cylinder(x, z, y);
    auto const y_cyl = transform_coordinate_cartesian_to_cylinder(y, z, y);
    auto const z_cyl = transform_coordinate_cartesian_to_cylinder(z, z, y);
    auto const x_ref = Vector3d{{1.0, -Utils::pi() / 2.0, 0.0}};
    auto const y_ref = Vector3d{{1.0, 0.0, 0.0}};
    auto const z_ref = Vector3d{{0.0, z_cyl[1], 1.0}};
    for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_SMALL(x_cyl[i] - x_ref[i], eps);
      BOOST_CHECK_SMALL(y_cyl[i] - y_ref[i], eps);
      BOOST_CHECK_SMALL(z_cyl[i] - z_ref[i], eps);
    }
  }
  // check transformation with orientation for another angle
  {
    auto const u = vec_rotate(z, Utils::pi() / 3.0, x);
    auto const v = vec_rotate(z, Utils::pi() / 3.0, y);
    auto const u_cyl = transform_coordinate_cartesian_to_cylinder(u, z, y);
    auto const v_cyl = transform_coordinate_cartesian_to_cylinder(v, z, y);
    auto const u_ref = Vector3d{{1.0, Utils::pi() * (1. / 3. - 1. / 2.), 0.0}};
    auto const v_ref = Vector3d{{1.0, Utils::pi() / 3.0, 0.0}};
    for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_SMALL(u_cyl[i] - u_ref[i], eps);
      BOOST_CHECK_SMALL(v_cyl[i] - v_ref[i], eps);
    }
  }
  // check transformation of random vectors
  {
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> rng(2);
    auto const r_uniform = [&rng]() {
      return static_cast<double>(rng() - rng.min()) / (rng.max() - rng.min());
    };
    for (int trial = 0; trial < 100; ++trial) {
      Vector3d const v1{r_uniform(), r_uniform(), r_uniform()};
      Vector3d const v2{r_uniform(), r_uniform(), r_uniform()};
      auto const a = Utils::vector_product(v1, v2) / v1.norm() / v2.norm();
      auto const v1_v1 = transform_coordinate_cartesian_to_cylinder(v1, a, v1);
      auto const v2_v1 = transform_coordinate_cartesian_to_cylinder(v2, a, v1);
      auto const v1_v2 = transform_coordinate_cartesian_to_cylinder(v1, a, v2);
      Vector3d const v1_v1_ref{v1.norm(), 0.0, 0.0};
      Vector3d const v2_v1_ref{v2.norm(), Utils::angle_between(v1, v2), 0.0};
      Vector3d const v1_v2_ref{v1.norm(), -Utils::angle_between(v1, v2), 0.0};
      for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_SMALL(v1_v1[i] - v1_v1_ref[i], eps);
        BOOST_CHECK_SMALL(v2_v1[i] - v2_v1_ref[i], eps);
        BOOST_CHECK_SMALL(v1_v2[i] - v1_v2_ref[i], eps);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(cylinder_to_cartesian_test) {
  constexpr auto eps = 1e-14;
  auto const cyl = Vector3d{{1.0, Utils::pi() / 4, 2.0}};
  auto const pos = transform_coordinate_cylinder_to_cartesian(cyl);
  BOOST_CHECK_SMALL(pos[0] - std::sqrt(2) / 2, eps);
  BOOST_CHECK_SMALL(pos[1] - std::sqrt(2) / 2, eps);
  BOOST_CHECK_SMALL(pos[2] - cyl[2], eps);
}

BOOST_AUTO_TEST_CASE(cylinder_to_cartesian_with_axis_and_orientation_test) {
  constexpr auto eps = 2e-14;
  Vector3d const cylinder_coord{{1.2, 3.123, 42.0}};
  auto const e_x = Vector3d{{1., 0., 0.}};
  auto const e_y = Vector3d{{0., 1., 0.}};
  auto const e_z = Vector3d{{0., 0., 1.}};

  auto const transformed_x =
      transform_coordinate_cylinder_to_cartesian(cylinder_coord, e_x, -e_z);
  auto const transformed_y =
      transform_coordinate_cylinder_to_cartesian(cylinder_coord, e_y, e_x);
  auto const transformed_z =
      transform_coordinate_cylinder_to_cartesian(cylinder_coord, e_z, e_x);
  // We transform from cylinder zu cartesian and have to rotate back. See test
  // cartesian_to_cylinder_test.
  auto const expected_x = vec_rotate(
      e_y, Utils::pi() / 2.0,
      transform_coordinate_cylinder_to_cartesian(cylinder_coord, e_z, e_x));
  auto const expected_y = vec_rotate(
      e_x, -Utils::pi() / 2.0,
      transform_coordinate_cylinder_to_cartesian(cylinder_coord, e_z, e_x));
  // x = r * cos(phi); y = r * sin(phi); z = z
  auto const expected_z = Vector3d{
      {cylinder_coord[0] * std::cos(cylinder_coord[1]),
       cylinder_coord[0] * std::sin(cylinder_coord[1]), cylinder_coord[2]}};
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK_SMALL(transformed_x[i] - expected_x[i], eps);
    BOOST_CHECK_SMALL(transformed_y[i] - expected_y[i], eps);
    BOOST_CHECK_SMALL(transformed_z[i] - expected_z[i], eps);
  }
}

BOOST_AUTO_TEST_CASE(cylinder_to_cartesian_with_axis_with_phi_2_test) {
  constexpr auto eps = 1e-14;
  // tilted orthogonal basis
  auto const y = (Vector3d{{0, 1, -1}}).normalize();
  auto const z = (Vector3d{{1, 1, 1}}).normalize();
  auto const x = Utils::vector_product(y, z);

  // check transformation with orientation
  {
    auto const x_cyl = transform_coordinate_cartesian_to_cylinder(x, z, y);
    auto const y_cyl = transform_coordinate_cartesian_to_cylinder(y, z, y);
    auto const z_cyl = transform_coordinate_cartesian_to_cylinder(z, z, y);
    auto const x_cart = transform_coordinate_cylinder_to_cartesian(x_cyl, z, y);
    auto const y_cart = transform_coordinate_cylinder_to_cartesian(y_cyl, z, y);
    auto const z_cart = transform_coordinate_cylinder_to_cartesian(z_cyl, z, y);
    for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_SMALL(x_cart[i] - x[i], eps);
      BOOST_CHECK_SMALL(y_cart[i] - y[i], eps);
      BOOST_CHECK_SMALL(z_cart[i] - z[i], eps);
    }
  }
  // check transformation with orientation for another angle
  {
    auto const u = vec_rotate(z, Utils::pi() / 3.0, x);
    auto const v = vec_rotate(z, Utils::pi() / 3.0, y);
    auto const u_cyl = transform_coordinate_cartesian_to_cylinder(u, z, y);
    auto const v_cyl = transform_coordinate_cartesian_to_cylinder(v, z, y);
    auto const u_cart = transform_coordinate_cylinder_to_cartesian(u_cyl, z, y);
    auto const v_cart = transform_coordinate_cylinder_to_cartesian(v_cyl, z, y);
    for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_SMALL(u_cart[i] - u[i], eps);
      BOOST_CHECK_SMALL(v_cart[i] - v[i], eps);
    }
  }
  // check transformation of random vectors
  {
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> rng(2);
    auto const r_uniform = [&rng]() {
      return static_cast<double>(rng() - rng.min()) / (rng.max() - rng.min());
    };
    for (int trial = 0; trial < 100; ++trial) {
      Vector3d const v1{r_uniform(), r_uniform(), r_uniform()};
      Vector3d const v2{r_uniform(), r_uniform(), r_uniform()};
      auto const a = Utils::vector_product(v1, v2) / v1.norm() / v2.norm();
      auto const v3 = transform_coordinate_cartesian_to_cylinder(v2, a, v1);
      auto const v4 = transform_coordinate_cylinder_to_cartesian(v3, a, v1);
      for (int i = 0; i < 3; ++i) {
        BOOST_CHECK_SMALL(v4[i] - v2[i], eps);
      }
    }
  }
}
