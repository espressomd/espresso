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

#include "reaction_methods/SingleReaction.hpp"

#include "reaction_methods/utils.hpp"

#include <cmath>
#include <unordered_map>

namespace ReactionMethods {

double factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i) {
  auto value = 1.;
  if (nu_i) {
    if (nu_i > 0) {
      for (int i = 1; i <= nu_i; i++) {
        value *= static_cast<double>(Ni0 + i);
      }
      value = 1. / value;
    } else {
      for (int i = 0; i < -nu_i; i++) {
        value *= static_cast<double>(Ni0 - i);
      }
    }
  }
  return value;
}

double calculate_factorial_expression(
    SingleReaction const &reaction,
    std::unordered_map<int, int> const &particle_numbers) {
  auto value = 1.;
  // factorial contribution of reactants
  for (int i = 0; i < reaction.reactant_types.size(); i++) {
    auto const nu_i = -1 * reaction.reactant_coefficients[i];
    auto const N_i0 = particle_numbers.at(reaction.reactant_types[i]);
    value *= factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N_i0, nu_i);
  }
  // factorial contribution of products
  for (int i = 0; i < reaction.product_types.size(); i++) {
    auto const nu_i = reaction.product_coefficients[i];
    auto const N_i0 = particle_numbers.at(reaction.product_types[i]);
    value *= factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N_i0, nu_i);
  }
  return value;
}

double calculate_factorial_expression_cpH(
    SingleReaction const &reaction,
    std::unordered_map<int, int> const &particle_numbers) {
  auto value = 1.;
  // factorial contribution of reactants
  {
    auto const nu_i = -1 * reaction.reactant_coefficients[0];
    auto const N_i0 = particle_numbers.at(reaction.reactant_types[0]);
    value *= factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N_i0, nu_i);
  }
  // factorial contribution of products
  {
    auto const nu_i = reaction.product_coefficients[0];
    auto const N_i0 = particle_numbers.at(reaction.product_types[0]);
    value *= factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(N_i0, nu_i);
  }
  return value;
}

} // namespace ReactionMethods
