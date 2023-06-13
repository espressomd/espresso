/*
 * Copyright (C) 2020-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE ReactionMethods utility functions test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "reaction_methods/utils.hpp"

#include <cmath>
#include <limits>

BOOST_AUTO_TEST_CASE(factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i_test) {
  using namespace ReactionMethods;
  auto constexpr tol = 8. * 100. * std::numeric_limits<double>::epsilon();

  auto const reaction_ensemble_combinations = [](int N, int nu) {
    return (N + nu < 0) ? 0. : std::tgamma(N + 1) / std::tgamma(N + nu + 1);
  };

  for (int N0 = 0; N0 < 6; ++N0) {
    for (int nu = -4; nu <= 4; ++nu) {
      auto const val = factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N0, nu);
      auto const ref = reaction_ensemble_combinations(N0, nu);
      BOOST_CHECK_CLOSE(val, ref, 10. * tol);
    }
  }
}
