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
#ifndef REACTION_METHODS_CONSTANT_PH_ENSEMBLE_HPP
#define REACTION_METHODS_CONSTANT_PH_ENSEMBLE_HPP

#include "reaction_methods/ExclusionRadius.hpp"
#include "reaction_methods/ReactionAlgorithm.hpp"

#include <map>
#include <memory>
#include <utility>

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
  ConstantpHEnsemble(boost::mpi::communicator const &comm, int seed, double kT,
                     std::shared_ptr<ExclusionRadius> exclusion,
                     double constant_pH)
      : ReactionAlgorithm(comm, seed, kT, std::move(exclusion)),
        m_constant_pH(constant_pH) {}
  double m_constant_pH;
};

} // namespace ReactionMethods
#endif
