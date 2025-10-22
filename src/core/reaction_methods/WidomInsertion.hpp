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
#ifndef REACTION_METHODS_WIDOM_INSERTION_HPP
#define REACTION_METHODS_WIDOM_INSERTION_HPP

#include "ReactionAlgorithm.hpp"

#include <unordered_map>
#include <utility>

namespace ReactionMethods {

/** Widom insertion method */
class WidomInsertion : public ReactionAlgorithm {
public:
  WidomInsertion(
      boost::mpi::communicator const &comm, int seed, double kT,
      double exclusion_range,
      const std::unordered_map<int, double> &exclusion_radius_per_type)
      : ReactionAlgorithm(comm, seed, kT, exclusion_range,
                          exclusion_radius_per_type) {}

  double calculate_particle_insertion_potential_energy(int reaction_id) {
    auto &reaction = *reactions[reaction_id];
    if (not all_reactant_particles_exist(reaction)) {
      throw std::runtime_error("Trying to remove some non-existing particles "
                               "from the system via the inverse Widom scheme.");
    }

    // make reaction attempt and immediately reverse it
    setup_bookkeeping_of_empty_pids();
    auto const E_pot_old = calculate_potential_energy();
    make_reaction_attempt(reaction, make_new_system_state());
    auto const E_pot_new = calculate_potential_energy();
    restore_old_system_state();

    return E_pot_new - E_pot_old;
  }
};

} // namespace ReactionMethods
#endif
