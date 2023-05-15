/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

#include "CylindricalLBVelocityProfile.hpp"

#include "grid_based_algorithms/lb_interface.hpp"

#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <algorithm>
#include <functional>
#include <vector>

namespace Observables {

std::vector<double> CylindricalLBVelocityProfile::operator()() const {
  Utils::CylindricalHistogram<double, 3> histogram(n_bins(), limits());
  for (auto const &p : sampling_positions) {
    auto const velocity =
        LB::get_interpolated_velocity(p) * LB::get_lattice_speed();
    auto const pos_shifted = p - transform_params->center();
    auto const pos_cyl = Utils::transform_coordinate_cartesian_to_cylinder(
        pos_shifted, transform_params->axis(), transform_params->orientation());
    histogram.update(pos_cyl,
                     Utils::transform_vector_cartesian_to_cylinder(
                         velocity, transform_params->axis(), pos_shifted));
  }
  auto hist_data = histogram.get_histogram();
  auto const tot_count = histogram.get_tot_count();
  std::transform(hist_data.begin(), hist_data.end(), tot_count.begin(),
                 hist_data.begin(), std::divides<double>());
  return hist_data;
}

} // namespace Observables
