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

#include "reaction_methods/ReactionEnsemble.hpp"

#include "utils.hpp"

#include <cmath>
#include <map>

namespace ReactionMethods {

/**
 * Calculates the expression in the acceptance probability in the reaction
 * ensemble
 */
double ReactionEnsemble::calculate_acceptance_probability(
    SingleReaction const &current_reaction, double E_pot_old, double E_pot_new,
    std::map<int, int> const &old_particle_numbers, int dummy_old_state_index,
    int dummy_new_state_index,
    bool dummy_only_make_configuration_changing_move) const {
  const double factorial_expr =
      calculate_factorial_expression(current_reaction, old_particle_numbers);

  const double beta = 1.0 / temperature;
  // calculate Boltzmann factor
  return std::pow(volume, current_reaction.nu_bar) * current_reaction.gamma *
         factorial_expr * exp(-beta * (E_pot_new - E_pot_old));
}

} // namespace ReactionMethods
