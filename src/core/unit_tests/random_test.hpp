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

/* Helper functions to compute random numbers covariance in a single pass */

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <boost/test/unit_test.hpp>
#include <tuple>

#include <utils/Vector.hpp>

#include "random.hpp"

/** Draw a large sample of 3D vectors from a PRNG and compute the following
 *  statistics: pooled mean, pooled variance, covariance matrix between axes.
 */
template <typename NoiseFunction>
std::tuple<double, double, Utils::Matrix<double, 3, 3>>
noise_stats(NoiseFunction noise_function, size_t sample_size) {
  namespace ba = boost::accumulators;
  namespace bt = boost::accumulators::tag;
  ba::accumulator_set<double, ba::stats<bt::mean, bt::variance(ba::lazy)>>
      acc_all, acc_x, acc_y, acc_z;
  ba::accumulator_set<double, ba::stats<bt::covariance<double, bt::covariate1>>>
      acc_xy, acc_xz, acc_yz;
  for (int i = 0; i < sample_size; ++i) {
    auto const noise = noise_function(i);
    acc_x(noise[0]);
    acc_y(noise[1]);
    acc_z(noise[2]);
    acc_xy(noise[0], ba::covariate1 = noise[1]);
    acc_xz(noise[0], ba::covariate1 = noise[2]);
    acc_yz(noise[1], ba::covariate1 = noise[2]);
    acc_all(noise[0]);
    acc_all(noise[1]);
    acc_all(noise[2]);
  }
  auto const mean = ba::mean(acc_all);
  auto const variance = ba::variance(acc_all);
  auto const cov_xx = ba::variance(acc_x);
  auto const cov_yy = ba::variance(acc_y);
  auto const cov_zz = ba::variance(acc_z);
  auto const cov_xy = ba::covariance(acc_xy);
  auto const cov_xz = ba::covariance(acc_xz);
  auto const cov_yz = ba::covariance(acc_yz);
  return std::make_tuple(mean, variance,
                         Utils::Matrix<double, 3, 3>{{cov_xx, cov_xy, cov_xz},
                                                     {cov_xy, cov_yy, cov_yz},
                                                     {cov_xz, cov_yz, cov_zz}});
}

void noise_check_stats(double expected_variance, double mean, double variance,
                       Utils::Matrix<double, 3, 3> const &cov) {
  auto const x = 0, y = 1, z = 2;
  // check total mean and variance
  BOOST_CHECK_SMALL(std::abs(mean), 1e-4);
  BOOST_CHECK_CLOSE(variance, expected_variance, 2e-2);
  // check variance per axis
  BOOST_CHECK_CLOSE(cov[x][x], expected_variance, 5e-2);
  BOOST_CHECK_CLOSE(cov[y][y], expected_variance, 5e-2);
  BOOST_CHECK_CLOSE(cov[z][z], expected_variance, 5e-2);
}

void noise_check_correlation(Utils::Matrix<double, 3, 3> const &cov) {
  auto const x = 0, y = 1, z = 2;
  // check the 3 axes are not correlated (Pearson correlation coefficient)
  auto const corrcoeff_xy = cov[x][y] / sqrt(cov[x][x] * cov[y][y]);
  auto const corrcoeff_xz = cov[x][z] / sqrt(cov[x][x] * cov[z][z]);
  auto const corrcoeff_yz = cov[y][z] / sqrt(cov[y][y] * cov[z][z]);
  BOOST_CHECK_SMALL(std::abs(corrcoeff_xy), 1e-3);
  BOOST_CHECK_SMALL(std::abs(corrcoeff_xz), 1e-3);
  BOOST_CHECK_SMALL(std::abs(corrcoeff_yz), 1e-3);
}
