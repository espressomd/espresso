/*
 * Copyright (C) 2016-2026 The ESPResSo project
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

#include "system/System.hpp"
#include "utils_histogram.hpp"

#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <vector>

namespace Observables {

std::vector<double> CylindricalLBVelocityProfile::operator()(
    boost::mpi::communicator const &comm) const {
  auto &system = System::get_system();
  auto &lb = system.lb;
  lb.ghost_communication_vel();

  if (lb_sanity_checks.mismatch(*system.box_geo, lb)) {
    calculate_sampling_positions(*system.box_geo, lb);
  }

  auto velocities = lb.get_coupling_interpolated_velocities(sampling_positions);
  auto pos_shifted_it = sampling_positions_cart.begin();
  std::vector<Utils::Vector3d> local_velocities{};
  local_velocities.reserve(velocities.size());

  for (auto const &vel : velocities) {
    auto const vel_cyl = Utils::transform_vector_cartesian_to_cylinder(
        vel, transform_params->axis(), *pos_shifted_it);
    local_velocities.emplace_back(vel_cyl);
    ++pos_shifted_it;
  }

  auto const [global_positions, global_velocities] =
      detail::gather(comm, sampling_positions_cyl, local_velocities);

  if (comm.rank() != 0) {
    return {};
  }

  Utils::CylindricalHistogram<double, 3> histogram(n_bins(), limits());
  detail::accumulate(histogram, global_positions, global_velocities);
  return detail::normalize_by_bin_size(histogram);
}

} // namespace Observables
