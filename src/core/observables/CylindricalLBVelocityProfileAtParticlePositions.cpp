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
#include "CylindricalLBVelocityProfileAtParticlePositions.hpp"

#include "BoxGeometry.hpp"
#include "system/System.hpp"
#include "utils_histogram.hpp"

#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <utility>
#include <vector>

namespace Observables {
std::vector<double> CylindricalLBVelocityProfileAtParticlePositions::evaluate(
    boost::mpi::communicator const &comm,
    ParticleReferenceRange const &local_particles,
    const ParticleObservables::traits<Particle> &traits) const {
  using pos_type = decltype(traits.position(std::declval<Particle>()));
  using vel_type = Utils::Vector3d;

  std::vector<pos_type> local_folded_positions{};
  std::vector<vel_type> local_velocities{};
  local_folded_positions.reserve(local_particles.size());
  local_velocities.reserve(local_particles.size());

  auto const &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto const &lb = system.lb;
  auto const vel_conv = lb.get_lattice_speed();

  for (auto const &p : local_particles) {
    auto const pos = box_geo.folded_position(traits.position(p));
    auto const pos_shifted = pos - transform_params->center();
    auto const vel = *lb.get_interpolated_velocity(pos);
    auto const pos_cyl = Utils::transform_coordinate_cartesian_to_cylinder(
        pos_shifted, transform_params->axis(), transform_params->orientation());
    auto const vel_cyl = Utils::transform_vector_cartesian_to_cylinder(
        vel * vel_conv, transform_params->axis(), pos_shifted);
    local_folded_positions.emplace_back(pos_cyl);
    local_velocities.emplace_back(vel_cyl);
  }

  auto const [global_folded_positions, global_velocities] =
      detail::gather(comm, local_folded_positions, local_velocities);

  if (comm.rank() != 0) {
    return {};
  }

  Utils::CylindricalHistogram<double, 3> histogram(n_bins(), limits());
  detail::accumulate(histogram, global_folded_positions, global_velocities);
  return detail::normalize_by_bin_size(histogram);
}

} // namespace Observables
