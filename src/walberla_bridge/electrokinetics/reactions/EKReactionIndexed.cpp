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

#include "EKReactionIndexed.hpp"
#include "EKReactant.hpp"
#include "EKReactionBase.hpp"
#include "LatticeWalberla.hpp"

#include "generated_kernels/ReactionKernelIndexed_1.h"
#include "generated_kernels/ReactionKernelIndexed_2.h"
#include "generated_kernels/ReactionKernelIndexed_3.h"
#include "generated_kernels/ReactionKernelIndexed_4.h"
#include "generated_kernels/ReactionKernelIndexed_5.h"

#include <memory>

#include "utils.hpp"

#include <field/AddToStorage.h>

namespace walberla {
namespace detail {
template <typename FloatType>
auto get_kernel(
    const std::vector<std::shared_ptr<EKReactant<FloatType>>> &reactants,
    const FloatType coefficient, const BlockDataID &indexFieldID) {
  switch (reactants.size()) {
  case 1: {
    const auto &reactant = reactants[0];
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactant);

    auto kernel = std::make_shared<pystencils::ReactionKernelIndexed_1>(
        indexFieldID, density_id_0, order_0, coefficient, stoech_coeff_0);

    return kernel->getSweep();
  }
  case 2: {
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactants[0]);
    const auto [density_id_1, order_1, stoech_coeff_1] =
        detail::get_reaction_details<FloatType>(reactants[1]);

    auto kernel = std::make_shared<pystencils::ReactionKernelIndexed_2>(
        indexFieldID, density_id_0, density_id_1, order_0, order_1, coefficient,
        stoech_coeff_0, stoech_coeff_1);

    return kernel->getSweep();
  }
  case 3: {
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactants[0]);
    const auto [density_id_1, order_1, stoech_coeff_1] =
        detail::get_reaction_details<FloatType>(reactants[1]);
    const auto [density_id_2, order_2, stoech_coeff_2] =
        detail::get_reaction_details<FloatType>(reactants[2]);

    auto kernel = std::make_shared<pystencils::ReactionKernelIndexed_3>(
        indexFieldID, density_id_0, density_id_1, density_id_2, order_0,
        order_1, order_2, coefficient, stoech_coeff_0, stoech_coeff_1,
        stoech_coeff_2);

    return kernel->getSweep();
  }
  case 4: {
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactants[0]);
    const auto [density_id_1, order_1, stoech_coeff_1] =
        detail::get_reaction_details<FloatType>(reactants[1]);
    const auto [density_id_2, order_2, stoech_coeff_2] =
        detail::get_reaction_details<FloatType>(reactants[2]);
    const auto [density_id_3, order_3, stoech_coeff_3] =
        detail::get_reaction_details<FloatType>(reactants[3]);

    auto kernel = std::make_shared<pystencils::ReactionKernelIndexed_4>(
        indexFieldID, density_id_0, density_id_1, density_id_2, density_id_3,
        order_0, order_1, order_2, order_3, coefficient, stoech_coeff_0,
        stoech_coeff_1, stoech_coeff_2, stoech_coeff_3);

    return kernel->getSweep();
  }
  case 5: {
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactants[0]);
    const auto [density_id_1, order_1, stoech_coeff_1] =
        detail::get_reaction_details<FloatType>(reactants[1]);
    const auto [density_id_2, order_2, stoech_coeff_2] =
        detail::get_reaction_details<FloatType>(reactants[2]);
    const auto [density_id_3, order_3, stoech_coeff_3] =
        detail::get_reaction_details<FloatType>(reactants[3]);
    const auto [density_id_4, order_4, stoech_coeff_4] =
        detail::get_reaction_details<FloatType>(reactants[4]);

    auto kernel = std::make_shared<pystencils::ReactionKernelIndexed_5>(
        indexFieldID, density_id_0, density_id_1, density_id_2, density_id_3,
        density_id_4, order_0, order_1, order_2, order_3, order_4, coefficient,
        stoech_coeff_0, stoech_coeff_1, stoech_coeff_2, stoech_coeff_3,
        stoech_coeff_4);

    return kernel->getSweep();
  }
  default:
    throw std::runtime_error("reactions of this size are not implemented!");
  }
}
} // namespace detail

template <typename FloatType>
EKReactionIndexed<FloatType>::EKReactionIndexed(
    std::shared_ptr<LatticeWalberla> lattice,
    std::vector<std::shared_ptr<EKReactant<FloatType>>> reactants,
    FloatType coefficient)
    : EKReactionBase<FloatType>(lattice, reactants, coefficient) {
  m_flagfield_id = field::addFlagFieldToStorage<FlagField>(
      get_lattice()->get_blocks(), "flag field reaction",
      get_lattice()->get_ghost_layers());

  using IndexVectors = pystencils::ReactionKernelIndexed_1::IndexVectors;

  auto createIdxVector = [](IBlock *const, StructuredBlockStorage *const) {
    return new IndexVectors();
  };
  m_indexvector_id = get_lattice()
                         ->get_blocks()
                         ->template addStructuredBlockData<IndexVectors>(
                             createIdxVector, "IndexField");
}

template <typename FloatType>
void EKReactionIndexed<FloatType>::perform_reaction() const {

  auto kernel =
      detail::get_kernel(get_reactants(), get_coefficient(), m_indexvector_id);

  for (auto &block : *get_lattice()->get_blocks()) {
    kernel(&block);
  }
}
} // namespace walberla