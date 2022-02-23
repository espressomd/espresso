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

#include "EKReactionImplBulk.hpp"
#include "EKReactant.hpp"

#include "generated_kernels/electrokinetics/reactions/ReactionKernel_1.h"
#include "generated_kernels/electrokinetics/reactions/ReactionKernel_2.h"
#include "generated_kernels/electrokinetics/reactions/ReactionKernel_3.h"
#include "generated_kernels/electrokinetics/reactions/ReactionKernel_4.h"
#include "generated_kernels/electrokinetics/reactions/ReactionKernel_5.h"

#include <domain_decomposition/BlockDataID.h>
#include <memory>

namespace walberla {

namespace detail {
template <typename FloatType>
auto get_reaction_details(
    const std::shared_ptr<EKReactant<FloatType>> &reactant) {
  const auto order = reactant->get_order();
  const auto stoech_coeff = reactant->get_stoech_coeff();
  const auto density_id =
      walberla::BlockDataID(reactant->get_species()->get_density_id());

  return std::make_tuple(density_id, order, stoech_coeff);
}

template <typename FloatType>
auto get_kernel(
    const std::vector<std::shared_ptr<EKReactant<FloatType>>> &reactants,
    const FloatType coefficient) {
  switch (reactants.size()) {
  case 1: {
    const auto &reactant = reactants[0];
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactant);

    auto kernel = std::make_shared<pystencils::ReactionKernel_1>(
        density_id_0, order_0, coefficient, stoech_coeff_0);

    return kernel->getSweep();
  }
  case 2: {
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactants[0]);
    const auto [density_id_1, order_1, stoech_coeff_1] =
        detail::get_reaction_details<FloatType>(reactants[1]);

    auto kernel = std::make_shared<pystencils::ReactionKernel_2>(
        density_id_0, density_id_1, order_0, order_1, coefficient,
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

    auto kernel = std::make_shared<pystencils::ReactionKernel_3>(
        density_id_0, density_id_1, density_id_2, order_0, order_1, order_2,
        coefficient, stoech_coeff_0, stoech_coeff_1, stoech_coeff_2);

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

    auto kernel = std::make_shared<pystencils::ReactionKernel_4>(
        density_id_0, density_id_1, density_id_2, density_id_3, order_0,
        order_1, order_2, order_3, coefficient, stoech_coeff_0, stoech_coeff_1,
        stoech_coeff_2, stoech_coeff_3);

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

    auto kernel = std::make_shared<pystencils::ReactionKernel_5>(
        density_id_0, density_id_1, density_id_2, density_id_3, density_id_4,
        order_0, order_1, order_2, order_3, order_4, coefficient,
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
void EKReactionImplBulk<FloatType>::perform_reaction() const {
  // TODO: if my understanding is correct:
  //  the kernels need to either run in the ghost layers and do the
  //  synchronization before or not run and do a synchronization afterwards.
  //  The better solution is probably the latter one. Not sure why it fails
  //  atm.

  auto kernel = detail::get_kernel(get_reactants(), get_coefficient());
  for (auto &block : *get_lattice()->get_blocks()) {
    kernel(&block);
  }
}
} // namespace walberla
