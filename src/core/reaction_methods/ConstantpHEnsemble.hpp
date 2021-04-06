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
#ifndef REACTION_METHODS_CONSTANT_PH_ENSEMBLE_HPP
#define REACTION_METHODS_CONSTANT_PH_ENSEMBLE_HPP

#include "reaction_methods/ReactionAlgorithm.hpp"

#include <map>

namespace ReactionMethods {

/**
 * Constant-pH Ensemble, for derivation see @cite reed92a.
 * For the constant pH reactions you need to provide the deprotonation and
 * afterwards the corresponding protonation reaction (in this order). If you
 * want to deal with multiple reactions do it multiple times. Note that there is
 * a difference in the usecase of the constant pH reactions and the above
 * reaction ensemble. For the constant pH simulation directily the
 * **apparent equilibrium constant which carries a unit** needs to be provided
 * -- this is equivalent to the gamma of the reaction ensemble above, where the
 * dimensionless reaction constant needs to be provided. Again: For the
 * constant-pH algorithm not the dimensionless reaction constant needs to be
 * provided here, but the apparent reaction constant.
 */
class ConstantpHEnsemble : public ReactionAlgorithm {
public:
  ConstantpHEnsemble(int seed) : ReactionAlgorithm(seed) {}
  double m_constant_pH = -10;
  int do_reaction(int reaction_steps) override;

protected:
  double calculate_acceptance_probability(
      SingleReaction const &current_reaction, double E_pot_old,
      double E_pot_new, std::map<int, int> const &dummy_old_particle_numbers,
      int dummy_old_state_index, int dummy_new_state_index,
      bool dummy_only_make_configuration_changing_move) const override;

private:
  int get_random_valid_p_id();
};

} // namespace ReactionMethods
#endif
