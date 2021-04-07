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

#include "reaction_methods/utils.hpp"

#include <cmath>
#include <map>

namespace ReactionMethods {

double factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i) {
  double value = 1.0;
  if (nu_i) {
    if (nu_i > 0) {
      for (int i = 1; i <= nu_i; i++) {
        value /= Ni0 + i;
      }
    } else {
      auto const abs_nu_i = std::abs(nu_i);
      for (int i = 0; i < abs_nu_i; i++) {
        value *= Ni0 - i;
      }
    }
  }
  return value;
}

double
calculate_factorial_expression(SingleReaction const &current_reaction,
                               std::map<int, int> const &old_particle_numbers) {
  double factorial_expr = 1.0;
  // factorial contribution of reactants
  for (int i = 0; i < current_reaction.reactant_types.size(); i++) {
    int nu_i = -1 * current_reaction.reactant_coefficients[i];
    int N_i0 = old_particle_numbers.at(current_reaction.reactant_types[i]);
    factorial_expr *= factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(
        N_i0, nu_i); // zeta = 1 (see @cite smith94c)
                     // since we only perform one reaction
                     // at one call of the function
  }
  // factorial contribution of products
  for (int i = 0; i < current_reaction.product_types.size(); i++) {
    int nu_i = current_reaction.product_coefficients[i];
    int N_i0 = old_particle_numbers.at(current_reaction.product_types[i]);
    factorial_expr *= factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(
        N_i0, nu_i); // zeta = 1 (see @cite smith94c)
                     // since we only perform one reaction
                     // at one call of the function
  }
  return factorial_expr;
}

} // namespace ReactionMethods
