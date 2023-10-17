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
#ifndef OBSERVABLES_FORCEDENSITYPROFILE_HPP
#define OBSERVABLES_FORCEDENSITYPROFILE_HPP

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "PidProfileObservable.hpp"
#include "system/System.hpp"
#include "utils_histogram.hpp"

#include <utils/Histogram.hpp>

#include <boost/mpi/collectives/gather.hpp>
#include <boost/serialization/vector.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace Observables {

class ForceDensityProfile : public PidProfileObservable {
public:
  using PidProfileObservable::PidProfileObservable;
  std::vector<std::size_t> shape() const override {
    auto const b = n_bins();
    return {b[0], b[1], b[2], 3};
  }

  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    using pos_type = decltype(traits.position(std::declval<Particle>()));
    using force_type = decltype(traits.force(std::declval<Particle>()));
    auto const &box_geo = *System::get_system().box_geo;

    std::vector<pos_type> local_folded_positions{};
    local_folded_positions.reserve(local_particles.size());
    std::vector<force_type> local_forces{};
    local_forces.reserve(local_particles.size());

    for (auto const &p : local_particles) {
      local_folded_positions.emplace_back(
          box_geo.folded_position(traits.position(p)));
      local_forces.emplace_back(traits.force(p));
    }

    auto const [global_folded_positions, global_forces] =
        detail::gather(comm, local_folded_positions, local_forces);

    if (comm.rank() != 0) {
      return {};
    }

    Utils::Histogram<double, 3> histogram(n_bins(), limits());
    detail::accumulate(histogram, global_folded_positions, global_forces);
    histogram.normalize();
    return histogram.get_histogram();
  }
};

} // Namespace Observables

#endif
