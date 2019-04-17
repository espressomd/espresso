/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define BOOST_TEST_MODULE AutoParameter test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "field_coupling/couplings/Charge.hpp"
#include "field_coupling/couplings/Direct.hpp"
#include "field_coupling/couplings/Mass.hpp"
#include "field_coupling/couplings/Scaled.hpp"
#include "field_coupling/couplings/Viscous.hpp"

using namespace FieldCoupling::Coupling;

BOOST_AUTO_TEST_CASE(charge) {
  static_assert(Charge::is_linear, "");

  struct {
    struct {
      const double q = 1.23;
    } p;
  } p;

  BOOST_CHECK((p.p.q * 2.0) == Charge()(p, 2.0));
}

BOOST_AUTO_TEST_CASE(mass) {
  static_assert(Mass::is_linear, "");

  struct {
    struct {
      const double mass = 1.23;
    } p;
  } p;

  BOOST_CHECK((p.p.mass * 2.0) == Mass()(p, 2.0));
}

BOOST_AUTO_TEST_CASE(direct) {
  static_assert(Direct::is_linear, "");

  BOOST_CHECK(5 == Direct()(0, 5));
}

BOOST_AUTO_TEST_CASE(scaled) {
  static_assert(Scaled::is_linear, "");

  auto const scales = std::unordered_map<int, double>{{0, 1.23}, {2, 3.45}};
  auto const default_val = 9.8;

  auto const scaled_coupling = Scaled(scales, default_val);

  BOOST_CHECK(scales == scaled_coupling.particle_scales());
  BOOST_CHECK(default_val == scaled_coupling.default_scale());

  {
    struct Particle {
      Particle(int id) : m_id(id) {}

      int identity() const { return m_id; };

      const int m_id;
    };

    BOOST_CHECK_CLOSE(1.23 * 2., scaled_coupling(Particle(0), 2.), 1e-14);
    BOOST_CHECK_CLOSE(default_val * 3., scaled_coupling(Particle(1), 3.),
                      1e-14);
    BOOST_CHECK_CLOSE(3.45 * 4., scaled_coupling(Particle(2), 4.), 1e-14);
    BOOST_CHECK_CLOSE(default_val * 5., scaled_coupling(Particle(3), 5.),
                      1e-14);
  }
}

BOOST_AUTO_TEST_CASE(viscous) {
  static_assert(Viscous::is_linear, "");

  auto const gamma = 3.14159;

  const auto viscous_coupling = Viscous(gamma);

  BOOST_CHECK(gamma == viscous_coupling.gamma());

  {
    struct {
      struct {
        const Utils::Vector3d v = {1., 2., 3.};
      } m;
    } p;

    auto const u = Utils::Vector3d{4., 5., 6.};

    BOOST_CHECK((-gamma * (p.m.v - u)) == viscous_coupling(p, u));
  }
}
