/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#define BOOST_TEST_MODULE Histogram test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/Histogram.hpp"
#include "utils/constants.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

BOOST_AUTO_TEST_CASE(histogram) {
  std::array<std::size_t, 2> n_bins{{10, 10}};
  std::array<std::pair<double, double>, 2> limits{
      {std::make_pair(1.0, 20.0), std::make_pair(5.0, 10.0)}};
  constexpr std::size_t n_dims_data = 2;
  auto hist = Utils::Histogram<double, n_dims_data, 2>(n_bins, limits);
  // Check getters.
  BOOST_CHECK(hist.get_limits() == limits);
  BOOST_CHECK(hist.get_n_bins() == n_bins);
  BOOST_CHECK((hist.get_bin_sizes() ==
               std::array<double, 2>{{19.0 / 10.0, 5.0 / 10.0}}));
  // Check that histogram is initialized to zero.
  BOOST_CHECK(hist.get_histogram() ==
              std::vector<double>(n_dims_data * n_bins[0] * n_bins[1], 0.0));
  // Check that histogram still empty if data is out of bounds.
  hist.update(std::vector<double>{{1.0, 4.0}});
  BOOST_CHECK(hist.get_histogram() ==
              std::vector<double>(n_dims_data * n_bins[0] * n_bins[1], 0.0));
  // Check if putting in data at the first bin is set correctly.
  hist.update(std::vector<double>{{limits[0].first, limits[1].first}});
  BOOST_CHECK((hist.get_histogram())[0] == 1.0);
  BOOST_CHECK((hist.get_histogram())[1] == 1.0);
  // Check if weights are correctly set.
  hist.update(std::vector<double>{{limits[0].first, limits[1].first}},
              std::vector<double>{{10.0, 10.0}});
  BOOST_CHECK((hist.get_histogram())[0] == 11.0);
  BOOST_CHECK((hist.get_histogram())[1] == 11.0);
  // Check exceptions
  BOOST_CHECK_THROW(hist.update(std::vector<double>{{1.0, 5.0, 3.0}}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(hist.update(std::vector<double>{{0.0, 0.0}},
                                std::vector<double>{{0.0, 0.0, 0.0}}),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(cylindrical_histogram) {
  constexpr auto pi = Utils::pi<double>();
  std::array<std::size_t, 3> n_bins{{10, 10, 10}};
  std::array<std::pair<double, double>, 3> limits{{std::make_pair(0.0, 2.0),
                                                   std::make_pair(0.0, 2 * pi),
                                                   std::make_pair(0.0, 10.0)}};
  constexpr std::size_t n_dims_data = 3;
  auto hist = Utils::CylindricalHistogram<double, n_dims_data>(n_bins, limits);
  // Check getters.
  BOOST_CHECK(hist.get_limits() == limits);
  BOOST_CHECK(hist.get_n_bins() == n_bins);
  BOOST_CHECK((hist.get_bin_sizes() ==
               std::array<double, 3>{{2.0 / 10.0, 2 * pi / 10.0, 1.0}}));
  // Check that histogram is initialized to zero.
  BOOST_CHECK(hist.get_histogram() ==
              std::vector<double>(n_dims_data * 1000, 0.0));
  // Check that histogram still empty if data is out of bounds.
  hist.update(std::vector<double>{{1.0, 3 * pi, 1.0}});
  BOOST_CHECK(hist.get_histogram() ==
              std::vector<double>(n_dims_data * 1000, 0.0));
  // Check if putting in data at the first bin is set correctly.
  hist.update(
      std::vector<double>{{limits[0].first, limits[1].first, limits[2].first}});
  BOOST_CHECK((hist.get_histogram())[0] == 1.0);
  BOOST_CHECK((hist.get_histogram())[1] == 1.0);
  BOOST_CHECK((hist.get_histogram())[2] == 1.0);
  // Check if weights are correctly set.
  hist.update(
      std::vector<double>{{limits[0].first, limits[1].first, limits[2].first}},
      std::vector<double>{{10.0, 10.0, 10.0}});
  BOOST_CHECK((hist.get_histogram())[0] == 11.0);
  BOOST_CHECK((hist.get_histogram())[1] == 11.0);
  BOOST_CHECK((hist.get_histogram())[2] == 11.0);
  // Check exceptions
  BOOST_CHECK_THROW(hist.update(std::vector<double>{{1.0, pi}}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(hist.update(std::vector<double>{{0.0, 0.0, 0.0}},
                                std::vector<double>{{0.0, 0.0}}),
                    std::invalid_argument);
}
