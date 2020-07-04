/*
 * Copyright (C) 2020 The ESPResSo project
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
#include "reaction_ensemble.hpp"

#include <limits>
#include <vector>

#define BOOST_TEST_MODULE ReactionEnsemble test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(find_minimum_non_negative_value_test) {
  using namespace ReactionEnsemble;

  BOOST_CHECK_EQUAL(find_minimum_non_negative_value({1, 2, 3}), 1);
  BOOST_CHECK_EQUAL(find_minimum_non_negative_value({3, 2, 1}), 1);
  BOOST_CHECK_EQUAL(find_minimum_non_negative_value({-1, 2, 3}), 2);
  BOOST_CHECK_EQUAL(find_minimum_non_negative_value({-1, -2, -3}), -3);
  BOOST_CHECK_THROW(find_minimum_non_negative_value({}), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(average_list_of_allowed_entries_test) {
  using namespace ReactionEnsemble;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  BOOST_CHECK_CLOSE(average_list_of_allowed_entries(std::vector<long>{1, 2}),
                    1.5, tol);
  BOOST_CHECK_CLOSE(
      average_list_of_allowed_entries(std::vector<double>{1, 2, -2}), 1.5, tol);
  BOOST_CHECK_CLOSE(
      average_list_of_allowed_entries(std::vector<double>{1.5, -3.}), 1.5, tol);
  BOOST_CHECK_CLOSE(
      average_list_of_allowed_entries(std::vector<int>{-1, -2, -3}), 0.0, tol);
  BOOST_CHECK_CLOSE(average_list_of_allowed_entries(std::vector<double>{}), 0.0,
                    tol);
}

BOOST_AUTO_TEST_CASE(factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i_test) {
  using namespace ReactionEnsemble;
  constexpr double tol = 100 * std::numeric_limits<double>::epsilon();

  BOOST_CHECK_CLOSE(factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(2, 0), 1.0,
                    tol);
  BOOST_CHECK_CLOSE(factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(2, 0), 1.0,
                    tol);
  BOOST_CHECK_CLOSE(factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(1, 2),
                    1.0 / 6.0, tol);
  BOOST_CHECK_CLOSE(factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(1, -2),
                    0.0, tol);
  BOOST_CHECK_CLOSE(factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(2, -2),
                    2.0, tol);
}
