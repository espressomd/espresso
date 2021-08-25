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

#define BOOST_TEST_MODULE sampling test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Histogram.hpp>
#include <utils/constants.hpp>
#include <utils/sampling.hpp>

#include <array>
#include <cstddef>
#include <utility>

BOOST_AUTO_TEST_CASE(get_cylindrical_sampling_positions_test) {
  auto const min_r = 0.0;
  auto const max_r = 5.0;
  auto const min_phi = -Utils::pi();
  auto const max_phi = Utils::pi();
  auto const min_z = 0.0;
  auto const max_z = 10.0;
  auto const n_r_bins = 10;
  auto const n_phi_bins = 10;
  auto const n_z_bins = 10;
  auto const sampling_density = 2.;
  auto const sampling_positions = Utils::get_cylindrical_sampling_positions(
      std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
      std::make_pair(min_z, max_z), n_r_bins, n_phi_bins, n_z_bins,
      sampling_density);
  std::array<std::pair<double, double>, 3> limits{
      {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
       std::make_pair(min_z, max_z)}};
  std::array<std::size_t, 3> n_bins{{static_cast<std::size_t>(n_r_bins),
                                     static_cast<std::size_t>(n_phi_bins),
                                     static_cast<std::size_t>(n_z_bins)}};
  Utils::CylindricalHistogram<double, 3> histogram(n_bins, limits);
  for (auto const &p : sampling_positions) {
    histogram.update(p);
  }
  auto const tot_count = histogram.get_tot_count();
  std::array<std::size_t, 3> const dimensions{{n_r_bins, n_phi_bins, n_z_bins}};
  std::array<std::size_t, 3> index{};
  for (auto const &c : tot_count) {
    BOOST_CHECK(c > 0);
  }
}
