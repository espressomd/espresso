/*
 * Copyright (C) 2016-2019 The ESPResSo project
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

#include <algorithm>

#include "CylindricalLBVelocityProfile.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

namespace Observables {

std::vector<double> CylindricalLBVelocityProfile::operator()() const {
  std::array<size_t, 3> n_bins{{n_r_bins, n_phi_bins, n_z_bins}};
  std::array<std::pair<double, double>, 3> limits{
      {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
       std::make_pair(min_z, max_z)}};
  Utils::CylindricalHistogram<double, 3> histogram(n_bins, 3, limits);
  for (auto const &p : sampling_positions) {
    auto const velocity = lb_lbfluid_get_interpolated_velocity(p) *
                          lb_lbfluid_get_lattice_speed();
    auto const pos_shifted = p - center;
    auto const pos_cyl =
        Utils::transform_coordinate_cartesian_to_cylinder(pos_shifted, axis);
    histogram.update(pos_cyl, Utils::transform_vector_cartesian_to_cylinder(
                                  velocity, axis, pos_shifted));
  }
  auto hist_data = histogram.get_histogram();
  auto const tot_count = histogram.get_tot_count();
  std::transform(hist_data.begin(), hist_data.end(), tot_count.begin(),
                 hist_data.begin(), std::divides<double>());
  return hist_data;
}

} // namespace Observables
