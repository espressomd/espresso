/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#define BOOST_TEST_MODULE observables_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <particle_observables/observable.hpp>

#include "mock.hpp"

#include <type_traits>
#include <vector>

namespace Testing {
template <class T> struct strip_args {
  template <class... Args> decltype(auto) operator()(Args...) const {
    return T{}();
  }
};
} // namespace Testing

BOOST_AUTO_TEST_CASE(product_) {
  using Testing::strip_args;

  auto prod = ParticleObservables::Product<
      strip_args<std::integral_constant<int, 2>>,
      strip_args<std::integral_constant<int, 3>>>{};

  BOOST_CHECK_EQUAL((prod.template operator()<int, int>(0)), 2 * 3);
}

BOOST_AUTO_TEST_CASE(obs) {
  using namespace ParticleObservables;
  Testing::Particle p;
  BOOST_CHECK_EQUAL(Momentum{}(p), Mass{}(p)*Velocity{}(p));
  std::vector<Testing::Particle> parts{p, Testing::Particle{}};
  parts[1].mass = 5.;
  {
    auto const res = AverageMomentum{}(parts).first;
    BOOST_CHECK(res == 0.5 * (Momentum{}(parts[0]) + Momentum{}(parts[1])));
  }
  {
    auto const res = CenterOfMassPosition{}(parts).first;
    BOOST_CHECK(res == (Mass{}(parts[0]) * Position{}(parts[0]) +
                        Mass{}(parts[1]) * Position{}(parts[1])) /
                           (Mass{}(parts[0]) + Mass{}(parts[1])));
  }
  {
    auto const res = CenterOfMassVelocity{}(parts).first;
    BOOST_CHECK(res == (Mass{}(parts[0]) * Velocity{}(parts[0]) +
                        Mass{}(parts[1]) * Velocity{}(parts[1])) /
                           (Mass{}(parts[0]) + Mass{}(parts[1])));
  }
  {
    auto const res = Positions{}(parts);
    BOOST_CHECK(res[0] == Position{}(parts[0]));
    BOOST_CHECK(res[1] == Position{}(parts[1]));
  }
  {
    auto const res = Velocities{}(parts);
    BOOST_CHECK(res[0] == Velocity{}(parts[0]));
    BOOST_CHECK(res[1] == Velocity{}(parts[1]));
  }
}
