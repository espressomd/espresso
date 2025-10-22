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
#ifndef OBSERVABLES_TotalForce_HPP
#define OBSERVABLES_TotalForce_HPP

#include "PidObservable.hpp"

#include <boost/mpi/collectives/reduce.hpp>

#include <cstddef>
#include <vector>

namespace Observables {
class TotalForce : public PidObservable {
public:
  using PidObservable::PidObservable;
  std::vector<std::size_t> shape() const override { return {3}; }

  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &) const override {
    Utils::Vector3d local_force{};
    for (auto const &p : local_particles) {
      if (p.get().is_virtual())
        continue;
      local_force += p.get().force();
    }

    decltype(local_force) global_force;
    boost::mpi::reduce(comm, local_force, global_force, std::plus<>(), 0);

    return global_force.as_vector();
  }
};
} // Namespace Observables
#endif
