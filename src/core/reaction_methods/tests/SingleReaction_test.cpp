/*
 * Copyright (C) 2021 The ESPResSo project
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

/* Unit tests for the SingleReaction. */

#define BOOST_TEST_MODULE SingleReaction test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "reaction_methods/SingleReaction.hpp"
#include "reaction_methods/utils.hpp"

#include <limits>
#include <map>
#include <memory>

// Check a simple chemical reaction, the Monte Carlo acceptance rate
// and the configurational move probability for a given system state.
BOOST_AUTO_TEST_CASE(SingleReaction_test) {
  using namespace ReactionMethods;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  // create a reaction A -> 3 B + 4 C
  int const type_A = 0;
  int const type_B = 1;
  int const type_C = 2;
  SingleReaction reaction(2., {type_A}, {1}, {type_B, type_C}, {3, 4});

  // check derived parameter nu_bar
  BOOST_CHECK_EQUAL(reaction.nu_bar, 6);

  // check acceptance rate
  for (int tried_moves = 1; tried_moves < 5; ++tried_moves) {
    for (int accepted_moves = 0; accepted_moves < 5; ++accepted_moves) {
      reaction.tried_moves = tried_moves;
      reaction.accepted_moves = accepted_moves;
      auto const ref_rate = static_cast<double>(accepted_moves) /
                            static_cast<double>(tried_moves);
      BOOST_CHECK_CLOSE(reaction.get_acceptance_rate(), ref_rate, tol);
    }
  }

  // check factorial expression
  constexpr auto g = factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        // system contains i x A, j x B, and k x C
        auto const p_numbers =
            std::map<int, int>{{type_A, i}, {type_B, j}, {type_C, k}};
        auto const val = calculate_factorial_expression(reaction, p_numbers);
        auto const ref = g(i, -1) * g(j, 3) * g(k, 4);
        BOOST_CHECK_CLOSE(val, ref, 5 * tol);
      }
    }
  }
}
