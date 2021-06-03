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
    std::map<int, int> const &old_particle_numbers, int, int, bool) const {
  auto const beta = 1.0 / temperature;
  double nu_H = 0;
  if (current_reaction.product_index_protons > 0)
    nu_H = current_reaction
               .product_coefficients[current_reaction.product_index_protons];
  else
    nu_H = -current_reaction
                .reactant_coefficients[-current_reaction.product_index_protons];

  auto const ln_bf = (E_pot_new - E_pot_old);
  const double factorial_expr = calculate_factorial_expression_cpH(
      current_reaction, old_particle_numbers);
  const double N_H = std::pow(10, -m_constant_pH) *
                     m_conversion_factor_from_mol_per_l_to_1_div_sigma_cubed *
                     volume;
  auto const bf = current_reaction.gamma *
                  std::pow(volume, current_reaction.nu_bar) *
                  std::pow(N_H, -nu_H) * factorial_expr * exp(-beta * ln_bf);
  return bf;
}

} // namespace ReactionMethods
