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
#include "LBVelocityProfile.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include <utils/Histogram.hpp>

namespace Observables {

std::vector<double> LBVelocityProfile::operator()() const {
  std::array<size_t, 3> n_bins{{n_x_bins, n_y_bins, n_z_bins}};
  std::array<std::pair<double, double>, 3> limits{
      {std::make_pair(min_x, max_x), std::make_pair(min_y, max_y),
       std::make_pair(min_z, max_z)}};
  Utils::Histogram<double, 3> histogram(n_bins, 3, limits);
  for (auto const &p : sampling_positions) {
    const auto v = lb_lbfluid_get_interpolated_velocity(p) *
                   lb_lbfluid_get_lattice_speed();
    histogram.update(p, v);
  }
  auto hist_tmp = histogram.get_histogram();
  auto const tot_count = histogram.get_tot_count();
  for (size_t ind = 0; ind < hist_tmp.size(); ++ind) {
    if (tot_count[ind] == 0 and not allow_empty_bins) {
      auto const error = "Decrease sampling delta(s), bin " +
                         std::to_string(ind) + " has no hit";
      throw std::runtime_error(error);
    }
    if (tot_count[ind] > 0) {
      hist_tmp[ind] /= tot_count[ind];
    }
  }
  return hist_tmp;
}

} // namespace Observables
