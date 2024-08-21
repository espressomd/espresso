/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "EKinWalberlaImpl.hpp"
#include "reactions/EKReactionImplBulk.hpp"
#include "reactions/EKReactionImplIndexed.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/ek_walberla_init.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBaseIndexed.hpp>

#include <utils/Vector.hpp>

#include <memory>

namespace walberla {

std::shared_ptr<EKinWalberlaBase>
new_ek_walberla(std::shared_ptr<LatticeWalberla> const &lattice,
                double diffusion, double kT, double valency,
                Utils::Vector3d ext_efield, double density, bool advection,
                bool friction_coupling, bool single_precision, bool thermalized,
                unsigned int seed) {
  if (single_precision) {
    return std::make_shared<EKinWalberlaImpl<13, float>>(
        lattice, diffusion, kT, valency, ext_efield, density, advection,
        friction_coupling, thermalized, seed);
  }

  return std::make_shared<EKinWalberlaImpl<13, double>>(
      lattice, diffusion, kT, valency, ext_efield, density, advection,
      friction_coupling, thermalized, seed);
}

std::shared_ptr<EKReactionBase>
new_ek_reaction_bulk(std::shared_ptr<LatticeWalberla> const &lattice,
                     typename EKReactionBase::reactants_type const &reactants,
                     double coefficient) {
  return std::make_shared<EKReactionImplBulk>(lattice, reactants, coefficient);
}

std::shared_ptr<EKReactionBaseIndexed> new_ek_reaction_indexed(
    std::shared_ptr<LatticeWalberla> const &lattice,
    typename EKReactionBase::reactants_type const &reactants,
    double coefficient) {
  return std::make_shared<EKReactionImplIndexed>(lattice, reactants,
                                                 coefficient);
}

} // namespace walberla
