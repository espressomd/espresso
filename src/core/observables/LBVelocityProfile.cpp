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
#include "utils_histogram.hpp"

#include <utils/Histogram.hpp>
#include <utils/Vector.hpp>

#include <boost/mpi/collectives/gather.hpp>
#include <boost/serialization/vector.hpp>

#include <stdexcept>
#include <vector>

namespace Observables {

std::vector<double>
LBVelocityProfile::operator()(boost::mpi::communicator const &comm) const {
  using vel_type = Utils::Vector3d;

  decltype(sampling_positions) local_positions{};
  std::vector<vel_type> local_velocities{};

  auto const vel_conv = LB::get_lattice_speed();

  for (auto const &pos : sampling_positions) {
    if (auto const vel = LB::get_interpolated_velocity(pos)) {
      local_positions.emplace_back(pos);
      local_velocities.emplace_back((*vel) * vel_conv);
    }
  }

  auto const world_size = comm.size();
  std::vector<decltype(sampling_positions)> global_positions{};
  std::vector<decltype(local_velocities)> global_velocities{};
  global_positions.reserve(world_size);
  global_velocities.reserve(world_size);
  boost::mpi::gather(comm, local_positions, global_positions, 0);
  boost::mpi::gather(comm, local_velocities, global_velocities, 0);

  if (comm.rank() != 0) {
    return {};
  }

  Utils::Histogram<double, 3> histogram(n_bins(), limits());
  accumulate(histogram, global_positions, global_velocities);
  try {
    return normalize_by_bin_size(histogram, allow_empty_bins);
  } catch (empty_bin_exception const &) {
    throw std::runtime_error(
        "Decrease sampling delta(s), some bins have no hit");
  }
}

} // namespace Observables
