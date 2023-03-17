/*
 * Copyright (C) 2022 The ESPResSo project
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

#define BOOST_TEST_MODULE rotation test
#define BOOST_TEST_DYN_LINK

#include "config/config.hpp"

#ifdef ROTATION

#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "rotation.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/quaternion.hpp>

#include <limits>
#include <stdexcept>
#include <tuple>

auto constexpr tol = 5. * 100. * std::numeric_limits<double>::epsilon();

namespace Testing {
std::tuple<Utils::Quaternion<double>, Utils::Vector3d>
setup_trivial_quat(int i, Utils::Vector3d const &v_in) {
  auto quat = Utils::Quaternion<double>{{0., 0., 0., 0.}};
  quat[i] = 1.;
  auto v_ref = v_in;
  if (i) {
    v_ref *= -1.;
    v_ref[i - 1] *= -1.;
  }
  return std::make_tuple(quat, v_ref);
}
} // namespace Testing

BOOST_AUTO_TEST_CASE(convert_vector_space_to_body_test) {
  auto const t_in = Utils::Vector3d{{1., 2., 3.}};
  for (int i : {0, 1, 2, 3}) {
    auto p = Particle();
    Utils::Vector3d t_ref;
    std::tie(p.quat(), t_ref) = Testing::setup_trivial_quat(i, t_in);
    auto const t_out = convert_vector_space_to_body(p, t_in);
    for (int j : {0, 1, 2}) {
      BOOST_CHECK_CLOSE(t_out[j], t_ref[j], tol);
    }
  }
}

BOOST_AUTO_TEST_CASE(convert_torque_to_body_frame_apply_fix_test) {
  auto const t_in = Utils::Vector3d{{1., 2., 3.}};
  {
    // test particle torque conversion
    for (int i : {0, 1, 2, 3}) {
      auto p = Particle();
      p.set_can_rotate_all_axes();
      Utils::Vector3d t_ref;
      std::tie(p.quat(), t_ref) = Testing::setup_trivial_quat(i, t_in);
      p.torque() = t_in;
      convert_torque_to_body_frame_apply_fix(p);
      auto const t_out = p.torque();
      for (int j : {0, 1, 2}) {
        BOOST_CHECK_CLOSE(t_out[j], t_ref[j], tol);
      }
    }
  }
  {
    // torque is set to zero for axes without rotation
    for (int j : {0, 1, 2}) {
      auto p = Particle();
      p.set_can_rotate_all_axes();
      p.set_can_rotate_around(j, false);
      p.quat() = Utils::Quaternion<double>::identity();
      p.torque() = t_in;
      auto t_ref = t_in;
      t_ref[j] = 0.;
      convert_torque_to_body_frame_apply_fix(p);
      BOOST_TEST(p.torque() == t_ref, boost::test_tools::per_element());
    }
  }
  {
    // torque is always zero for non-rotatable particles
    auto p = Particle();
    p.set_cannot_rotate_all_axes();
    p.quat() = Utils::Quaternion<double>::identity();
    p.torque() = t_in;
    convert_torque_to_body_frame_apply_fix(p);
    auto const t_ref = Utils::Vector3d{};
    BOOST_TEST(p.torque() == t_ref, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_CASE(rotate_particle_body_test) {
  auto p = Particle();
  p.quat() = {1., 2., 3., 4.};
  {
    // fixed particles are unaffected, quaternion is identical to original
    p.set_cannot_rotate_all_axes();
    auto const phi = 2.;
    auto const quat = local_rotate_particle_body(p, {0., 0., 1.}, phi);
    BOOST_TEST((quat == p.quat()));
  }
  {
    // edge case: null rotation throws an exception
    p.set_can_rotate_around(2, true);
    auto const phi = 2.;
    BOOST_CHECK_THROW(local_rotate_particle_body(p, {1., 1., 0.}, phi),
                      std::exception);
  }
  {
    // an angle of zero has no effect, quaternion is identical to original
    p.set_can_rotate_all_axes();
    auto const phi = 0.;
    auto const quat = local_rotate_particle_body(p, {1., 2., 3.}, phi);
    BOOST_TEST((quat == p.quat()));
  }
  {
    // an angle of pi around the z-axis flips the quaternion sequence
    p.set_can_rotate_all_axes();
    auto const phi = Utils::pi<double>();
    auto const quat = local_rotate_particle_body(p, {0., 0., 1.}, phi);
    auto const quat_ref = Utils::Vector4d{{-4., 3., -2., 1.}};
    for (int i : {0, 1, 2, 3}) {
      BOOST_CHECK_CLOSE(quat[i], quat_ref[i], tol);
    }
  }
  {
    // an angle of 2 pi around the z-axis flips the quaternion sign
    p.set_can_rotate_all_axes();
    auto const phi = 2. * Utils::pi<double>();
    auto const quat = local_rotate_particle_body(p, {0., 0., 1.}, phi);
    auto const quat_ref = Utils::Vector4d{{-1., -2., -3., -4.}};
    for (int i : {0, 1, 2, 3}) {
      BOOST_CHECK_CLOSE(quat[i], quat_ref[i], tol);
    }
  }
}

BOOST_AUTO_TEST_CASE(propagate_omega_quat_particle_test) {
  auto p = Particle();
  p.set_can_rotate_all_axes();
  {
    // test edge case: null quaternion and no rotation
    p.quat() = {0., 0., 0., 0.};
    p.omega() = {0., 0., 0.};
    propagate_omega_quat_particle(p, 0.01);
    auto const quat = p.quat();
    auto const quat_ref = Utils::Quaternion<double>::identity();
    for (int i : {0, 1, 2, 3}) {
      BOOST_CHECK_CLOSE(quat[i], quat_ref[i], tol);
    }
  }
  {
    // test trivial cases with parameters extremely close to the limit:
    // time step almost 1.0 and product of time step with omega almost 2.0
    auto const time_step = 0.99;
    for (int j : {0, 1, 2}) {
      p.quat() = Utils::Quaternion<double>::identity();
      p.omega() = {0., 0., 0.};
      p.omega()[j] = 2.;
      propagate_omega_quat_particle(p, time_step);
      auto const quat = p.quat();
      auto quat_ref = Utils::Quaternion<double>::identity();
      quat_ref[1 + j] = time_step;
      quat_ref[0] = std::sqrt(1. - time_step * time_step);
      for (int i : {0, 1, 2, 3}) {
        BOOST_CHECK_CLOSE(quat[i], quat_ref[i], tol);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(convert_operator_body_to_space_test) {
  auto constexpr sqrt_2_half = Utils::sqrt_2() / 2.0;
  // rotation around y-axis by pi/2
  Utils::Quaternion<double> const quat = {sqrt_2_half, 0.0, sqrt_2_half, 0.0};
  // rotation around z-axis by pi/4
  Utils::Matrix<double, 3, 3> const linear_transf_body = {
      {sqrt_2_half, -sqrt_2_half, 0.0},
      {sqrt_2_half, sqrt_2_half, 0.0},
      {0.0, 0.0, 1.0}};
  // rotation around x-axis by pi/4
  Utils::Matrix<double, 3, 3> const linear_transf_space_ref = {
      {1.0, 0.0, 0.0},
      {0.0, sqrt_2_half, -sqrt_2_half},
      {0.0, sqrt_2_half, sqrt_2_half}};

  auto const linear_transf_space =
      convert_body_to_space(quat, linear_transf_body);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (linear_transf_space_ref(i, j) == 0.0) {
        BOOST_CHECK_SMALL(
            std::abs(linear_transf_space(i, j) - linear_transf_space_ref(i, j)),
            tol);
      } else {
        BOOST_CHECK_CLOSE(linear_transf_space(i, j),
                          linear_transf_space_ref(i, j), tol);
      }
    }
  }
}

#ifdef DIPOLES
BOOST_AUTO_TEST_CASE(convert_dip_to_quat_test) {
  auto const quat_to_vector4d = [](Utils::Quaternion<double> const &quat) {
    return Utils::Vector4d{quat.data(), quat.data() + 4};
  };
  auto p = Particle();
  p.quat() = {1., 2., 3., 4.};
  {
    auto const dipm = 0.8;
    auto const pair = convert_dip_to_quat({0., 0., dipm});
    auto const quat = quat_to_vector4d(pair.first);
    auto const quat_ref = Utils::Vector4d{{1., 0., 0., 0.}};
    for (int i : {0, 1, 2, 3}) {
      BOOST_CHECK_CLOSE(quat[i], quat_ref[i], tol);
    }
    BOOST_CHECK_CLOSE(pair.second, dipm, tol);
  }
  {
    auto const dipm = 1.6;
    auto const pair = convert_dip_to_quat({dipm, 0., 0.});
    auto const quat = quat_to_vector4d(pair.first);
    auto const quat_ref = Utils::Vector4d{{0.5, -0.5, 0.5, -0.5}};
    for (int i : {0, 1, 2, 3}) {
      BOOST_CHECK_CLOSE(quat[i], quat_ref[i], tol);
    }
    BOOST_CHECK_CLOSE(pair.second, dipm, tol);
  }
}
#endif // DIPOLES
#else  // ROTATION
int main(int argc, char **argv) {}
#endif // ROTATION
