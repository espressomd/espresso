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

#include "utils.hpp"

#include <cmath>
#include <map>
#include <vector>

namespace ReactionMethods {

/**
 * Calculates the expression in the acceptance probability of the constant pH
 * method.
 */
double ConstantpHEnsemble::calculate_acceptance_probability(
    SingleReaction const &current_reaction, double E_pot_old, double E_pot_new,
    std::map<int, int> const &old_particle_numbers, int dummy_old_state_index,
    int dummy_new_state_index,
    bool dummy_only_make_configuration_changing_move) const {
  auto const beta = 1.0 / temperature;
  auto const pKa = -current_reaction.nu_bar * log10(current_reaction.gamma);
  auto const ln_bf = (E_pot_new - E_pot_old) - current_reaction.nu_bar / beta *
                                                   log(10) *
                                                   (m_constant_pH - pKa);
  const double factorial_expr = calculate_factorial_expression_cpH(
      current_reaction, old_particle_numbers);
  auto const bf = factorial_expr * exp(-beta * ln_bf);
  return bf;
}

} // namespace ReactionMethods
