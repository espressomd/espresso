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

#include "BoxGeometry.hpp"
#include "LeesEdwardsBC.hpp"
#include "Particle.hpp"
#include "lees_edwards.hpp"

#include <utils/Vector.hpp>

#include <boost/range/algorithm/equal.hpp>

#include <algorithm>
#include <limits>

using namespace LeesEdwards;
using Utils::Vector3d;

BOOST_AUTO_TEST_CASE(test_shear_direction) {
  LeesEdwardsBC le;
  le.shear_direction = 1;
  BoxGeometry box;
  box.set_lees_edwards_bc(le);
  BOOST_CHECK_SMALL((shear_direction(box) - Vector3d{{0, 1, 0}}).norm2(),
                    std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(test_update_offset) {
  Particle p;
  p.l.i = {2, 4, -1};
  p.l.lees_edwards_offset = 1.5;

  LeesEdwardsBC le;
  le.shear_direction = 1;
  le.shear_plane_normal = 2;
  le.shear_velocity = 3.5;
  BoxGeometry box;
  box.set_lees_edwards_bc(le);
  update_offset(p, box, 3.5);
  BOOST_CHECK_CLOSE(p.l.lees_edwards_offset,
                    1.5 - le.shear_velocity * 0.5 * 3.5 *
                              (p.l.i[le.shear_plane_normal]),
                    std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(test_push) {
  const double old_offset = -1.2;
  const double shear_l = 6;
  const double shear_normal_l = 4.5;
  const double dt = 0.8;
  const Vector3d old_pos{3, shear_normal_l * 1.1, 10};
  const Vector3d old_vel{-1.2, 2, 4.1};

  Particle p;

  p.r.p = old_pos;
  p.m.v = old_vel;

  p.l.i = {2, 4, -1};
  p.l.lees_edwards_offset = old_offset;

  LeesEdwardsBC le;
  le.shear_direction = 2;
  le.shear_plane_normal = 1;
  le.pos_offset = 2.5;
  le.shear_velocity = -3.1;

  BoxGeometry box;
  box.set_type(BoxType::LEES_EDWARDS);
  box.set_length({5, shear_normal_l, shear_l});
  box.set_lees_edwards_bc(le);

  push(p, box, dt);

  auto expected_pos = old_pos - shear_direction(box) * le.pos_offset;
  BOOST_CHECK_SMALL((p.r.p - expected_pos).norm2(),
                    std::numeric_limits<double>::epsilon());

  auto expected_vel = old_vel - shear_direction(box) * le.shear_velocity;
  BOOST_CHECK_SMALL((p.m.v - expected_vel).norm2(),
                    std::numeric_limits<double>::epsilon());
  BOOST_CHECK_CLOSE(p.l.lees_edwards_offset,
                    old_offset + le.pos_offset -
                        le.shear_velocity * 0.5 * dt *
                            p.l.i[le.shear_plane_normal],
                    std::numeric_limits<double>::epsilon());

  // Test transition in the other diretion
  p.r.p[le.shear_plane_normal] = -1;
  expected_pos = {old_pos[0], -1, expected_pos[2]};
  expected_vel = old_vel;
  push(p, box, dt);
  BOOST_CHECK_CLOSE(p.l.lees_edwards_offset,
                    old_offset -
                        le.shear_velocity * dt * p.l.i[le.shear_plane_normal],
                    std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(protocols) {
  auto off = Off();
  BOOST_CHECK_EQUAL(get_pos_offset(17.3, off), 0.0);
  BOOST_CHECK_EQUAL(get_shear_velocity(17.3, off), 0.0);

  const double x0 = -2.1;
  const double t0 = 1.2;
  const double v = 2.6;

  auto linear = LinearShear(x0, v, t0);
  BOOST_CHECK_CLOSE(get_pos_offset(3.3, linear), x0 + v * (3.3 - t0),
                    std::numeric_limits<double>::epsilon());
  BOOST_CHECK_EQUAL(get_shear_velocity(17.3, linear), v);

  const double a = 3.1;
  const double o = 2.1;
  auto osc = OscillatoryShear(a, o, t0);
  BOOST_CHECK_CLOSE(get_pos_offset(3.3, osc), a * sin(o * (3.3 - t0)),
                    std::numeric_limits<double>::epsilon());
  BOOST_CHECK_CLOSE(get_shear_velocity(3.3, osc), a * o * cos(o * (3.3 - t0)),
                    std::numeric_limits<double>::epsilon());
}
