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
#ifndef OBSERVABLES_CYLINDRICALFLUXDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALFLUXDENSITYPROFILE_HPP

#include "BoxGeometry.hpp"
#include "CylindricalPidProfileObservable.hpp"
#include "grid.hpp"
#include "utils_histogram.hpp"

#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace Observables {
class CylindricalFluxDensityProfile : public CylindricalPidProfileObservable {
public:
  using CylindricalPidProfileObservable::CylindricalPidProfileObservable;

  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    using pos_type = decltype(traits.position(std::declval<Particle>()));
    using vel_type = decltype(traits.velocity(std::declval<Particle>()));

    std::vector<pos_type> local_folded_positions{};
    local_folded_positions.reserve(local_particles.size());
    std::vector<vel_type> local_velocities{};
    local_velocities.reserve(local_particles.size());

    for (auto const &p : local_particles) {
      auto const pos = folded_position(traits.position(p), box_geo) -
                       transform_params->center();

      local_folded_positions.emplace_back(
          Utils::transform_coordinate_cartesian_to_cylinder(
              pos, transform_params->axis(), transform_params->orientation()));
      local_velocities.emplace_back(
          Utils::transform_vector_cartesian_to_cylinder(
              traits.velocity(p), transform_params->axis(), pos));
    }

    auto const world_size = comm.size();
    std::vector<decltype(local_folded_positions)> global_folded_positions{};
    std::vector<decltype(local_velocities)> global_velocities{};
    global_folded_positions.reserve(world_size);
    global_velocities.reserve(world_size);
    boost::mpi::gather(comm, local_folded_positions, global_folded_positions,
                       0);
    boost::mpi::gather(comm, local_velocities, global_velocities, 0);

    if (comm.rank() != 0) {
      return {};
    }

    Utils::CylindricalHistogram<double, 3> histogram(n_bins(), limits());
    accumulate(histogram, global_folded_positions, global_velocities);
    histogram.normalize();
    return histogram.get_histogram();
  }
  std::vector<std::size_t> shape() const override {
    auto const b = n_bins();
    return {b[0], b[1], b[2], 3};
  }
};

} // Namespace Observables

#endif
