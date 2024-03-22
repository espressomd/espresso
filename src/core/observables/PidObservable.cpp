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
#include "PidObservable.hpp"

#include "ParticleTraits.hpp"
#include "fetch_particles.hpp"

#include <boost/mpi/communicator.hpp>

#include <vector>

namespace Observables {
std::vector<double>
PidObservable::operator()(boost::mpi::communicator const &comm) const {
  auto const &local_particles = fetch_particles(ids());
  return this->evaluate(comm, local_particles,
                        ParticleObservables::traits<Particle>{});
}
} // namespace Observables
