/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include <cmath>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace ReactionMethods {

double WidomInsertion::calculate_particle_insertion_potential_energy(
    SingleReaction &current_reaction) {

  if (!all_reactant_particles_exist(current_reaction))
    throw std::runtime_error("Trying to remove some non-existing particles "
                             "from the system via the inverse Widom scheme.");

  auto const E_pot_old = calculate_current_potential_energy_of_system();

  // make reaction attempt
  std::vector<int> p_ids_created_particles;
  std::vector<StoredParticleProperty> hidden_particles_properties;
  std::vector<StoredParticleProperty> changed_particles_properties;

  // save p_id, charge and type of the reactant particle, only thing we
  // need to hide the particle and recover it
  auto constexpr number_of_saved_properties = 3;

  std::tie(changed_particles_properties, p_ids_created_particles,
           hidden_particles_properties) =
      make_reaction_attempt(current_reaction);

  auto const E_pot_new = calculate_current_potential_energy_of_system();
  // reverse reaction attempt
  // reverse reaction
  // 1) delete created product particles
  for (int p_ids_created_particle : p_ids_created_particles) {
    delete_particle(p_ids_created_particle);
  }
  // 2) restore previously hidden reactant particles
  restore_properties(hidden_particles_properties, number_of_saved_properties);
  // 3) restore previously changed reactant particles
  restore_properties(changed_particles_properties, number_of_saved_properties);

  // calculate the particle insertion potential energy
  auto const E_pot_insertion = E_pot_new - E_pot_old;

  return E_pot_insertion;
}

} // namespace ReactionMethods
