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

BOOST_AUTO_TEST_CASE(test_noise_statistics) {
  constexpr size_t const sample_size = 100'000;
  constexpr size_t const x = 0, y = 1, z = 2;
  constexpr double const tol = 1e-12;

  double value = 1;
  std::vector<double> means, variances;
  std::vector<std::vector<double>> covariance;
  std::vector<std::vector<double>> correlation;
  std::tie(means, variances, covariance, correlation) =
      noise_statistics(std::function<std::vector<VariantVectorXd>()>(
                           [&value]() -> std::vector<VariantVectorXd> {
                             value *= -1;
                             return {{Utils::Vector2d{value, -value}}};
                           }),
                       sample_size);
  // check pooled mean and variance
  BOOST_CHECK_SMALL(std::abs(means[0]), 100 * tol);
  BOOST_CHECK_CLOSE(variances[0], 1.0, tol);
  // check variance per axis
  BOOST_CHECK_CLOSE(covariance[x][x], 1.0, tol);
  BOOST_CHECK_CLOSE(covariance[y][y], 1.0, tol);
  BOOST_CHECK_CLOSE(covariance[x][y], -1.0, tol);
  BOOST_CHECK_EQUAL(covariance[x][y], covariance[y][x]);
  // check correlation
  BOOST_CHECK_CLOSE(correlation[x][x], 1.0, tol);
  BOOST_CHECK_CLOSE(correlation[y][y], 1.0, tol);
  BOOST_CHECK_CLOSE(correlation[x][y], -1.0, tol);
  BOOST_CHECK_EQUAL(correlation[x][y], correlation[y][x]);
}

BOOST_AUTO_TEST_CASE(test_noise_uniform_1d) {
  constexpr size_t const sample_size = 4'000'000;

  std::vector<double> means, variances;
  std::vector<std::vector<double>> covariance;
  std::vector<std::vector<double>> correlation;
  std::tie(means, variances, covariance, correlation) = noise_statistics(
      std::function<std::vector<VariantVectorXd>()>(
          [counter = 0]() mutable -> std::vector<VariantVectorXd> {
            return {{Random::noise_uniform<RNGSalt::NPTISOV, 1>(counter++, 0)}};
          }),
      sample_size);
  // check pooled mean and variance
  BOOST_CHECK_SMALL(std::abs(means[0]), 2e-4);
  BOOST_CHECK_CLOSE(variances[0] * 12.0, 1.0, 0.05);
}

BOOST_AUTO_TEST_CASE(test_noise_uniform_3d) {
  constexpr size_t const sample_size = 4'000'000;
  constexpr size_t const x = 0, y = 1, z = 2;

  std::vector<double> means, variances;
  std::vector<std::vector<double>> covariance;
  std::vector<std::vector<double>> correlation;
  std::tie(means, variances, covariance, correlation) = noise_statistics(
      std::function<std::vector<VariantVectorXd>()>(
          [counter = 0]() mutable -> std::vector<VariantVectorXd> {
            return {{Random::noise_uniform<RNGSalt::LANGEVIN>(counter++, 0)}};
          }),
      sample_size);
  // check pooled mean and variance
  BOOST_CHECK_SMALL(std::abs(means[0]), 2e-4);
  BOOST_CHECK_CLOSE(variances[0] * 12.0, 1.0, 0.05);
  // check variance per axis
  BOOST_CHECK_CLOSE(covariance[x][x] * 12.0, 1.0, 0.2);
  BOOST_CHECK_CLOSE(covariance[y][y] * 12.0, 1.0, 0.2);
  BOOST_CHECK_CLOSE(covariance[z][z] * 12.0, 1.0, 0.2);
  // check correlation
  BOOST_CHECK_SMALL(std::abs(correlation[x][y]), 2e-3);
  BOOST_CHECK_SMALL(std::abs(correlation[y][z]), 2e-3);
  BOOST_CHECK_SMALL(std::abs(correlation[z][x]), 2e-3);
}

BOOST_AUTO_TEST_CASE(test_noise_gaussian_4d) {
  constexpr size_t const sample_size = 5'000'000;
  constexpr size_t const x = 0, y = 1, z = 2, t = 3;

  std::vector<double> means, variances;
  std::vector<std::vector<double>> covariance;
  std::vector<std::vector<double>> correlation;
  std::tie(means, variances, covariance, correlation) = noise_statistics(
      std::function<std::vector<VariantVectorXd>()>(
          [counter = 0]() mutable -> std::vector<VariantVectorXd> {
            return {{Random::noise_gaussian<RNGSalt::BROWNIAN_WALK, 4>(
                counter++, 0)}};
          }),
      sample_size);
  // check pooled mean and variance
  BOOST_CHECK_SMALL(std::abs(means[0]), 2e-4);
  BOOST_CHECK_CLOSE(variances[0], 1.0, 0.05);
  // check variance per axis
  BOOST_CHECK_CLOSE(covariance[x][x], 1.0, 0.2);
  BOOST_CHECK_CLOSE(covariance[y][y], 1.0, 0.2);
  BOOST_CHECK_CLOSE(covariance[z][z], 1.0, 0.2);
  BOOST_CHECK_CLOSE(covariance[t][t], 1.0, 0.2);
  // check correlation
  BOOST_CHECK_SMALL(std::abs(correlation[x][y]), 2e-3);
  BOOST_CHECK_SMALL(std::abs(correlation[y][z]), 2e-3);
  BOOST_CHECK_SMALL(std::abs(correlation[z][x]), 2e-3);
  BOOST_CHECK_SMALL(std::abs(correlation[x][t]), 2e-3);
  BOOST_CHECK_SMALL(std::abs(correlation[y][t]), 2e-3);
  BOOST_CHECK_SMALL(std::abs(correlation[z][t]), 2e-3);
}
