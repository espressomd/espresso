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
#include <utility>
#include <vector>

namespace ReactionMethods {

std::pair<double, double>
WidomInsertion::measure_excess_chemical_potential(int reaction_id) {
  if (!all_reactant_particles_exist(reaction_id))
    throw std::runtime_error("Trying to remove some non-existing particles "
                             "from the system via the inverse Widom scheme.");

  SingleReaction &current_reaction = reactions[reaction_id];
  const double E_pot_old = calculate_current_potential_energy_of_system();

  // make reaction attempt
  std::vector<int> p_ids_created_particles;
  std::vector<StoredParticleProperty> hidden_particles_properties;
  std::vector<StoredParticleProperty> changed_particles_properties;
  const int number_of_saved_properties =
      3; // save p_id, charge and type of the reactant particle, only thing we
         // need to hide the particle and recover it
  make_reaction_attempt(current_reaction, changed_particles_properties,
                        p_ids_created_particles, hidden_particles_properties);
  const double E_pot_new = calculate_current_potential_energy_of_system();
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
  std::vector<double> exponential = {
      exp(-1.0 / temperature * (E_pot_new - E_pot_old))};
  current_reaction.accumulator_exponentials(exponential);

  // calculate mean excess chemical potential and standard error of the mean
  std::pair<double, double> result = std::make_pair(
      -temperature * log(current_reaction.accumulator_exponentials.mean()[0]),
      std::abs(-temperature /
               current_reaction.accumulator_exponentials.mean()[0] *
               current_reaction.accumulator_exponentials.std_error()[0]));
  return result;
}

} // namespace ReactionMethods
