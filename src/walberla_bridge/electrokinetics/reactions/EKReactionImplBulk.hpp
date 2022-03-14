/*
 * Copyright (C) 2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_WALBERLA_BRIDGE_EKREACTIONIMPLBULK_HPP
#define ESPRESSO_SRC_WALBERLA_BRIDGE_EKREACTIONIMPLBULK_HPP

#include "EKReactant.hpp"
#include "EKReactionBase.hpp"
#include "LatticeWalberla.hpp"

#include <memory>

namespace walberla {

template <typename FloatType>
class EKReactionImplBulk : public EKReactionBase<FloatType> {
private:
  using ReactionBase = EKReactionBase<FloatType>;

public:
  EKReactionImplBulk(
      std::shared_ptr<LatticeWalberla> lattice,
      std::vector<std::shared_ptr<EKReactant<FloatType>>> reactants,
      FloatType coefficient)
      : EKReactionBase<FloatType>(lattice, reactants, coefficient) {}

  using ReactionBase::get_coefficient;
  using ReactionBase::get_lattice;
  using ReactionBase::get_reactants;

  void perform_reaction() override;
};

// explicit template instantiation
template class EKReactionImplBulk<double>;
} // namespace walberla

#endif // ESPRESSO_SRC_WALBERLA_BRIDGE_EKREACTIONIMPLBULK_HPP
