/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace Observables {

std::vector<double> LBVelocityProfile::operator()() const {
  Utils::Histogram<double, 3> histogram(n_bins(), limits());
  for (auto const &p : sampling_positions) {
    const auto v = LB::get_interpolated_velocity(p) * LB::get_lattice_speed();
    histogram.update(p, v);
  }
  auto hist_tmp = histogram.get_histogram();
  auto const tot_count = histogram.get_tot_count();
  for (std::size_t ind = 0; ind < hist_tmp.size(); ++ind) {
    if (tot_count[ind] == 0 and not allow_empty_bins) {
      auto const error = "Decrease sampling delta(s), bin " +
                         std::to_string(ind) + " has no hit";
      throw std::runtime_error(error);
    }
    if (tot_count[ind] > 0) {
      hist_tmp[ind] /= static_cast<double>(tot_count[ind]);
    }
  }
  return hist_tmp;
}

} // namespace Observables
