/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

#include "random.hpp"
#include "random_test.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <tuple>
#include <vector>

BOOST_AUTO_TEST_CASE(test_noise_statistics) {
  constexpr std::size_t const sample_size = 60'000;
  constexpr std::size_t const x = 0, y = 1;
  constexpr double const tol = 1e-12;

  double value = 1;
  auto const [means, variances, covariance, correlation] = noise_statistics(
      [&value]() -> std::array<VariantVectorXd, 1> {
        value *= -1;
        return {{Utils::Vector2d{value, -value}}};
      },
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
  constexpr std::size_t const sample_size = 60'000;

  auto const [means, variances, covariance, correlation] = noise_statistics(
      [counter = 0]() mutable -> std::array<VariantVectorXd, 1> {
        return {{Random::noise_uniform<RNGSalt::NPTISOV, 1>(counter++, 0, 1)}};
      },
      sample_size);
  // check pooled mean and variance
  BOOST_CHECK_SMALL(std::abs(means[0]), 1e-3);
  BOOST_CHECK_CLOSE(variances[0] * 12.0, 1.0, 0.1);
}

BOOST_AUTO_TEST_CASE(test_noise_uniform_3d) {
  constexpr std::size_t const sample_size = 60'000;
  constexpr std::size_t const x = 0, y = 1, z = 2;

  auto const [means, variances, covariance, correlation] = noise_statistics(
      [counter = 0]() mutable -> std::array<VariantVectorXd, 1> {
        return {{Random::noise_uniform<RNGSalt::LANGEVIN>(counter++, 1, 0)}};
      },
      sample_size);
  // check pooled mean and variance
  BOOST_CHECK_SMALL(std::abs(means[0]), 1e-3);
  BOOST_CHECK_CLOSE(variances[0] * 12.0, 1.0, 0.1);
  // check variance per axis
  BOOST_CHECK_CLOSE(covariance[x][x] * 12.0, 1.0, 0.8);
  BOOST_CHECK_CLOSE(covariance[y][y] * 12.0, 1.0, 0.8);
  BOOST_CHECK_CLOSE(covariance[z][z] * 12.0, 1.0, 0.8);
  // check correlation
  BOOST_CHECK_SMALL(std::abs(correlation[x][y]), 1.1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[y][z]), 1.1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[z][x]), 1.1e-2);
}

BOOST_AUTO_TEST_CASE(test_noise_gaussian_4d) {
  constexpr std::size_t const sample_size = 100'000;
  constexpr std::size_t const x = 0, y = 1, z = 2, t = 3;

  auto const [means, variances, covariance, correlation] = noise_statistics(
      [counter = 0]() mutable -> std::array<VariantVectorXd, 1> {
        return {{Random::noise_gaussian<RNGSalt::BROWNIAN_WALK, 4>(counter++, 0,
                                                                   0)}};
      },
      sample_size);
  // check pooled mean and variance
  BOOST_CHECK_SMALL(std::abs(means[0]), 1e-3);
  BOOST_CHECK_CLOSE(variances[0], 1.0, 0.15);
  // check variance per axis
  BOOST_CHECK_CLOSE(covariance[x][x], 1.0, 0.6);
  BOOST_CHECK_CLOSE(covariance[y][y], 1.0, 0.6);
  BOOST_CHECK_CLOSE(covariance[z][z], 1.0, 0.6);
  BOOST_CHECK_CLOSE(covariance[t][t], 1.0, 0.6);
  // check correlation
  BOOST_CHECK_SMALL(std::abs(correlation[x][y]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[y][z]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[z][x]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[x][t]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[y][t]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[z][t]), 1e-2);
}

BOOST_AUTO_TEST_CASE(test_uncorrelated_consecutive_ids) {
  // setup: 2 particles with the same seed and consecutive ids
  // check thermostats with pid offset by 2 aren't cross-correlated with lag 2
  constexpr std::size_t const sample_size = 50'000;
  constexpr std::size_t const x = 0, y = 1, z = 2;
  constexpr std::size_t seed = 0;
  constexpr int pid = 1;
  constexpr int pid_offset = 2;

  auto const correlation = std::get<3>(noise_statistics(
      [counter = 0]() mutable -> std::array<VariantVectorXd, 3> {
        counter++;
        auto prng = Random::noise_uniform<RNGSalt::NPTISOV, 1>;
        return {{prng(counter, seed, pid, 0),
                 prng(counter, seed, pid + pid_offset, 0),
                 prng(counter + pid_offset, seed, pid, 0)}};
      },
      sample_size));
  // check correlation
  BOOST_CHECK_SMALL(std::abs(correlation[x][y]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[x][z]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[y][z]), 1e-2);
}

BOOST_AUTO_TEST_CASE(test_uncorrelated_consecutive_seeds) {
  // setup: 2 particles with the same id with 2 rngs with consecutive seeds
  // check thermostats with seed offset by 2 aren't cross-correlated with lag 2
  constexpr std::size_t const sample_size = 50'000;
  constexpr std::size_t const x = 0, y = 1, z = 2;
  constexpr int pid = 1;
  constexpr std::size_t seed = 0;
  constexpr std::size_t seed_offset = 2;

  auto const correlation = std::get<3>(noise_statistics(
      [counter = 0]() mutable -> std::array<VariantVectorXd, 3> {
        counter++;
        auto prng = Random::noise_uniform<RNGSalt::NPTISOV, 1>;
        return {{prng(counter, seed, pid, 0),
                 prng(counter, seed + seed_offset, pid, 0),
                 prng(counter + seed_offset, seed, pid, 0)}};
      },
      sample_size));
  // check correlation
  BOOST_CHECK_SMALL(std::abs(correlation[x][y]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[x][z]), 1e-2);
  BOOST_CHECK_SMALL(std::abs(correlation[y][z]), 1e-2);
}
