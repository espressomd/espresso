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

#include "system/System.hpp"
#include "utils_histogram.hpp"

#include <utils/Histogram.hpp>
#include <utils/Vector.hpp>

#include <stdexcept>
#include <vector>

namespace Observables {

std::vector<double>
LBVelocityProfile::operator()(boost::mpi::communicator const &comm) const {
  using vel_type = Utils::Vector3d;

  decltype(sampling_positions) local_positions{};
  std::vector<vel_type> local_velocities{};

  auto const &lb = System::get_system().lb;
  auto const vel_conv = lb.get_lattice_speed();

  for (auto const &pos : sampling_positions) {
    if (auto const vel = lb.get_interpolated_velocity(pos)) {
      local_positions.emplace_back(pos);
      local_velocities.emplace_back((*vel) * vel_conv);
    }
  }

  auto const [global_positions, global_velocities] =
      detail::gather(comm, local_positions, local_velocities);

  if (comm.rank() != 0) {
    return {};
  }

  Utils::Histogram<double, 3> histogram(n_bins(), limits());
  detail::accumulate(histogram, global_positions, global_velocities);
  try {
    return detail::normalize_by_bin_size(histogram, allow_empty_bins);
  } catch (detail::empty_bin_exception const &) {
    throw std::runtime_error(
        "Decrease sampling delta(s), some bins have no hit");
  }
}

} // namespace Observables
