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

#include "reaction_methods/ConstantpHEnsemble.hpp"

#include "particle_data.hpp"

#include <cmath>
#include <map>
#include <vector>

namespace ReactionMethods {

int ConstantpHEnsemble::get_random_valid_p_id() {
  int random_p_id = i_random(get_maximal_particle_id() + 1);
  // draw random p_ids till we draw a pid which exists
  while (not particle_exists(random_p_id))
    random_p_id = i_random(get_maximal_particle_id() + 1);
  return random_p_id;
}

/**
 *Performs a reaction in the constant pH ensemble
 */
int ConstantpHEnsemble::do_reaction(int reaction_steps) {

  for (int i = 0; i < reaction_steps; ++i) {
    // get a list of reactions where a randomly selected particle type occurs in
    // the reactant list. the selection probability of the particle types has to
    // be proportional to the number of occurrences of the number of particles
    // with
    // this type

    // for optimizations this list could be determined during the initialization
    std::vector<int> list_of_reaction_ids_with_given_reactant_type;
    while (list_of_reaction_ids_with_given_reactant_type
               .empty()) { // avoid selecting a (e.g. salt) particle which
                           // does not take part in a reaction
      int random_p_id = get_random_valid_p_id(); // only used to determine which
                                                 // reaction is attempted.
      auto part = get_particle_data(random_p_id);

      int type_of_random_p_id = part.p.type;

      // construct list of reactions with the above reactant type
      for (int reaction_i = 0; reaction_i < reactions.size(); reaction_i++) {
        SingleReaction &current_reaction = reactions[reaction_i];
        for (int reactant_i = 0; reactant_i < 1;
             reactant_i++) { // reactant_i < 1 since it is assumed in this place
                             // that the types A and HA occur in the first place
                             // only. These are the types that should be
                             // switched, H+ should not be switched
          if (current_reaction.reactant_types[reactant_i] ==
              type_of_random_p_id) {
            list_of_reaction_ids_with_given_reactant_type.push_back(reaction_i);
            break;
          }
        }
      }
    }
    // randomly select a reaction to be performed
    int reaction_id =
        list_of_reaction_ids_with_given_reactant_type[i_random(static_cast<int>(
            list_of_reaction_ids_with_given_reactant_type.size()))];
    generic_oneway_reaction(reaction_id);
  }
  return 0;
}

/**
 * Calculates the expression in the acceptance probability of the constant pH
 * method.
 */
double ConstantpHEnsemble::calculate_acceptance_probability(
    SingleReaction const &current_reaction, double E_pot_old, double E_pot_new,
    std::map<int, int> const &dummy_old_particle_numbers,
    int dummy_old_state_index, int dummy_new_state_index,
    bool dummy_only_make_configuration_changing_move) const {
  auto const beta = 1.0 / temperature;
  auto const pKa = -current_reaction.nu_bar * log10(current_reaction.gamma);
  auto const ln_bf = (E_pot_new - E_pot_old) - current_reaction.nu_bar / beta *
                                                   log(10) *
                                                   (m_constant_pH - pKa);
  auto const bf = exp(-beta * ln_bf);
  return bf;
}

} // namespace ReactionMethods
