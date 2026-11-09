/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
#include "cell_system/CellStructure.hpp"
#include "fetch_particles.hpp"
#include "system/System.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>

#include <functional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Observables {

/** Validate that every requested particle id corresponds to an existing
 *  particle. The documented precondition for an observable is that its ids
 *  refer to existing particles; otherwise @ref fetch_particles silently drops
 *  the missing ids, which downstream leads to an out-of-bounds read in the
 *  chain observables and silently-wrong results in the map observables.
 *
 *  The check is collective so it cannot deadlock: every rank participates in
 *  the same reductions and throws the same exception, which is surfaced to
 *  Python via the script interface's parallel_try_catch wrapper.
 */
static void check_particle_ids_exist(boost::mpi::communicator const &comm,
                                     std::vector<int> const &ids) {
  auto const &cell_structure = *System::get_system().cell_structure;
  auto local_count = 0;
  for (auto const pid : ids) {
    auto const *ptr = cell_structure.get_local_particle(pid);
    if (ptr != nullptr and not ptr->is_ghost()) {
      ++local_count;
    }
  }
  auto const global_count =
      boost::mpi::all_reduce(comm, local_count, std::plus<int>{});
  if (global_count != static_cast<int>(ids.size())) {
    // Identify the first missing id (collectively) for a helpful message.
    for (auto const pid : ids) {
      auto const *ptr = cell_structure.get_local_particle(pid);
      auto const present = (ptr != nullptr and not ptr->is_ghost()) ? 1 : 0;
      auto const total =
          boost::mpi::all_reduce(comm, present, std::plus<int>{});
      if (total == 0) {
        std::stringstream error_msg;
        error_msg << "Particle with id " << pid << " does not exist; "
                  << "the observable cannot be computed because it references "
                  << "particle ids that do not exist.";
        throw std::runtime_error(error_msg.str());
      }
    }
  }
}

std::vector<double>
PidObservable::operator()(boost::mpi::communicator const &comm) const {
  check_particle_ids_exist(comm, ids());
  auto const &local_particles = fetch_particles(ids());
  return this->evaluate(comm, local_particles,
                        ParticleObservables::traits<Particle>{});
}
} // namespace Observables
