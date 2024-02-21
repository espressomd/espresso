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

#pragma once

#include "generated_kernels/ReactionKernelBulk_all.h"

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactant.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp>

#include <blockforest/StructuredBlockForest.h>

namespace walberla {

class EKReactionImplBulk : public EKReactionBase {
public:
  ~EKReactionImplBulk() override = default;

  using EKReactionBase::EKReactionBase;
  using EKReactionBase::get_coefficient;
  using EKReactionBase::get_lattice;
  using EKReactionBase::get_reactants;

  void perform_reaction() override {
    // TODO: if my understanding is correct:
    // the kernels need to either run in the ghost layers and do the
    // synchronization before or not run and do a synchronization afterwards.
    // The better solution is probably the latter one. Not sure why it fails
    // atm.
    auto kernel = detail::ReactionKernelBulkSelector::get_kernel(
        get_reactants(), get_coefficient());
    for (auto &block : *get_lattice()->get_blocks()) {
      kernel(&block);
    }
  }
};

} // namespace walberla
