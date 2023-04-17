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
#define BOOST_TEST_MODULE algorithms_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <particle_observables/algorithms.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

using namespace ParticleObservables;

namespace Testing {
struct Identity {
  template <class Particle> auto operator()(Particle const &p) { return p; }
};
struct One {
  template <class Particle> auto operator()(Particle const &) { return 1; }
};
struct PlusOne {
  template <class Particle> auto operator()(Particle const &p) { return p + 1; }
};

template <typename T>
T average(std::vector<T> const &value, std::size_t denominator) {
  auto const sum = std::accumulate(value.begin(), value.end(), T{0});
  return static_cast<T>(sum) / static_cast<T>(denominator);
}

} // namespace Testing

BOOST_AUTO_TEST_CASE(algorithms_integer) {
  std::vector<int> const values{1, 2, 3, 4};
  {
    auto const res = WeightedAverage<Testing::Identity, Testing::One>()(values);
    BOOST_CHECK_EQUAL(res, Testing::average(values, values.size()));
  }
  {
    auto const res =
        WeightedAverage<Testing::Identity, Testing::PlusOne>()(values);
    BOOST_CHECK_EQUAL(res, (1 * 2 + 2 * 3 + 3 * 4 + 4 * 5) / 14);
  }
  {
    auto const res = WeightedSum<Testing::Identity, Testing::PlusOne>()(values);
    BOOST_CHECK_EQUAL(res, (1 * 2 + 2 * 3 + 3 * 4 + 4 * 5));
  }
  {
    auto const res = Average<Testing::Identity>()(values);
    BOOST_CHECK_EQUAL(res, Testing::average(values, values.size()));
  }
  {
    auto const res = Sum<Testing::Identity>{}(values);
    BOOST_CHECK_EQUAL(res, Testing::average(values, 1u));
  }
  {
    auto const res = Map<Testing::Identity>{}(values);
    BOOST_TEST(res == values);
  }
}

BOOST_AUTO_TEST_CASE(algorithms_double) {
  auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();
  std::vector<double> const values{1., 2., 3., 4.};
  {
    auto const res = WeightedAverage<Testing::Identity, Testing::One>()(values);
    BOOST_CHECK_CLOSE(res, Testing::average(values, values.size()), tol);
  }
  {
    auto const res = Average<Testing::Identity>()(values);
    BOOST_CHECK_EQUAL(res, Testing::average(values, values.size()));
  }
  {
    auto const res = Sum<Testing::Identity>{}(values);
    BOOST_CHECK_EQUAL(res, Testing::average(values, 1u));
  }
  {
    auto const res = Map<Testing::Identity>{}(values);
    BOOST_TEST(res == values);
  }
}
