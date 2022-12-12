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

#include "reaction_methods/WidomInsertion.hpp"

#include "energy.hpp"

#include <stdexcept>

namespace ReactionMethods {

double WidomInsertion::calculate_particle_insertion_potential_energy(
    SingleReaction &current_reaction) {

  if (!all_reactant_particles_exist(current_reaction))
    throw std::runtime_error("Trying to remove some non-existing particles "
                             "from the system via the inverse Widom scheme.");

  auto const E_pot_old = mpi_calculate_potential_energy();

  // Setup the list of empty pids for bookeeping
  setup_bookkeeping_of_empty_pids();

  // make reaction attempt and immediately reverse it
  auto const change_tracker = make_reaction_attempt(current_reaction);
  auto const E_pot_new = mpi_calculate_potential_energy();
  change_tracker.restore_original_state();

  // calculate the particle insertion potential energy
  auto const E_pot_insertion = E_pot_new - E_pot_old;

  return E_pot_insertion;
}

} // namespace ReactionMethods
