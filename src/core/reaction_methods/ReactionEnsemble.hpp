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
#ifndef REACTION_METHODS_REACTION_ENSEMBLE_HPP
#define REACTION_METHODS_REACTION_ENSEMBLE_HPP

#include "reaction_methods/ReactionAlgorithm.hpp"

#include <map>

namespace ReactionMethods {

/** Reaction ensemble method.
 *  Works for the reaction ensemble at constant volume and temperature. For the
 *  reaction ensemble at constant pressure, additionally employ a barostat!
 *  NOTE: a chemical reaction consists of a forward and backward reaction.
 *  Here both reactions have to be defined separately. The extent of the
 *  reaction is here chosen to be +1. If the reaction trial move for a
 *  dissociation of HA is accepted then there is one more dissociated ion
 *  pair H+ and A-. Implementation of @cite smith94c.
 */
class ReactionEnsemble : public ReactionAlgorithm {
public:
  ReactionEnsemble(int seed) : ReactionAlgorithm(seed) {}

protected:
  double calculate_acceptance_probability(
      SingleReaction const &current_reaction, double E_pot_old,
      double E_pot_new, std::map<int, int> const &old_particle_numbers,
      int dummy_old_state_index, int dummy_new_state_index,
      bool dummy_only_make_configuration_changing_move) const override;
};

} // namespace ReactionMethods
#endif
