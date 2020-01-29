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

/* Unit tests for random number generators. */

#define BOOST_TEST_MODULE PRNG test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <tuple>

#include <utils/Vector.hpp>

#include "random.hpp"
#include "random_test.hpp"

BOOST_AUTO_TEST_CASE(test_noise_uniform) {
  double mean, var;
  Utils::Matrix<double, 3, 3> cov;
  std::tie(mean, var, cov) = noise_stats(
      [](int i) -> Utils::Vector3d {
        return Random::v_noise<RNGSalt::LANGEVIN>(i, 0);
      },
      10'000'000);
  noise_check_stats(1.0 / 12.0, mean, var, cov);
  noise_check_correlation(cov);
}

BOOST_AUTO_TEST_CASE(test_noise_gaussian) {
  double mean, var;
  Utils::Matrix<double, 3, 3> cov;
  std::tie(mean, var, cov) = noise_stats(
      [](int i) -> Utils::Vector3d {
        return Random::v_noise_g<RNGSalt::BROWNIAN_WALK>(i, 0);
      },
      10'000'000);
  noise_check_stats(1.0, mean, var, cov);
  noise_check_correlation(cov);
}
