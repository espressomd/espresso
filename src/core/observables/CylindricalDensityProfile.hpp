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
#ifndef OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP

#include "CylindricalPidProfileObservable.hpp"

#include "BoxGeometry.hpp"
#include "system/System.hpp"
#include "utils_histogram.hpp"

#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

#include <cstddef>
#include <vector>

namespace Observables {
class CylindricalDensityProfile : public CylindricalPidProfileObservable {
public:
  using CylindricalPidProfileObservable::CylindricalPidProfileObservable;
  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    using pos_type = Utils::Vector3d;
    auto const &box_geo = *System::get_system().box_geo;

    std::vector<pos_type> local_folded_positions{};
    local_folded_positions.reserve(local_particles.size());

    for (auto const &p : local_particles) {
      auto const pos = box_geo.folded_position(traits.position(p));
      auto const pos_shifted = pos - transform_params->center();

      local_folded_positions.emplace_back(
          Utils::transform_coordinate_cartesian_to_cylinder(
              pos_shifted, transform_params->axis(),
              transform_params->orientation()));
    }

    auto const global_folded_positions =
        detail::gather(comm, local_folded_positions);

    if (comm.rank() != 0) {
      return {};
    }

    Utils::CylindricalHistogram<double, 1> histogram(n_bins(), limits());

    for (auto const &vec : global_folded_positions) {
      for (auto const &p : vec) {
        histogram.update(p);
      }
    }

    histogram.normalize();
    return histogram.get_histogram();
  }
};

} // Namespace Observables

#endif
