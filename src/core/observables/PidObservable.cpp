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
#include <boost/mpi/operations.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace Observables {

void check_particle_ids_exist(boost::mpi::communicator const &comm,
                              std::vector<int> const &ids) {
  auto const &cell_structure = *System::get_system().cell_structure;
  auto const is_local = [&cell_structure](int pid) {
    auto const *p = cell_structure.get_local_particle(pid);
    return p != nullptr and not p->is_ghost();
  };
  auto const local_count = std::ranges::count_if(ids, is_local);
  auto const global_count =
      boost::mpi::all_reduce(comm, local_count, std::plus<>{});
  if (global_count == std::ranges::ssize(ids)) {
    return;
  }
  std::vector<int> local_present(ids.size());
  std::ranges::transform(ids, local_present.begin(), [&is_local](int pid) {
    return is_local(pid) ? 1 : 0;
  });
  std::vector<int> global_present(ids.size());
  boost::mpi::all_reduce(comm, local_present.data(),
                         static_cast<int>(ids.size()), global_present.data(),
                         boost::mpi::maximum<int>{});
  auto const missing = std::ranges::find(global_present, 0);
  assert(missing != global_present.end());
  auto const pid = ids[static_cast<std::size_t>(
      std::distance(global_present.begin(), missing))];
  throw std::runtime_error("Particle with id " + std::to_string(pid) +
                           " does not exist");
}

std::vector<double>
PidObservable::operator()(boost::mpi::communicator const &comm) const {
  check_particle_ids_exist(comm, ids());
  auto const &local_particles = fetch_particles(ids());
  return this->evaluate(comm, local_particles,
                        ParticleObservables::traits<Particle>{});
}
} // namespace Observables
