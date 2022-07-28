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
#ifndef REACTION_METHODS_UTILS_HPP
#define REACTION_METHODS_UTILS_HPP

#include "reaction_methods/ReactionAlgorithm.hpp"

#include <map>

namespace ReactionMethods {

/**
 * Calculates the whole product of factorial expressions which occur in the
 * reaction ensemble acceptance probability.
 *
 * See @cite smith94c.
 */
double
calculate_factorial_expression(SingleReaction const &current_reaction,
                               std::map<int, int> const &old_particle_numbers);

/**
 * Calculates the factorial expression which occurs in the constant pH method
 * with symmetric proposal probability.
 *
 * See @cite landsgesell17b
 */
double calculate_factorial_expression_cpH(
    SingleReaction const &current_reaction,
    std::map<int, int> const &old_particle_numbers);

/**
 * Calculates the factorial expression which occurs in the reaction ensemble
 * acceptance probability
 */
double factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i);

} // namespace ReactionMethods
#endif
