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

#include "generated_kernels/ReactionKernelBulk_1_double_precision.h"
#include "generated_kernels/ReactionKernelBulk_1_single_precision.h"
#include "generated_kernels/ReactionKernelBulk_2_double_precision.h"
#include "generated_kernels/ReactionKernelBulk_2_single_precision.h"
#include "generated_kernels/ReactionKernelBulk_3_double_precision.h"
#include "generated_kernels/ReactionKernelBulk_3_single_precision.h"
#include "generated_kernels/ReactionKernelBulk_4_double_precision.h"
#include "generated_kernels/ReactionKernelBulk_4_single_precision.h"
#include "generated_kernels/ReactionKernelBulk_5_double_precision.h"
#include "generated_kernels/ReactionKernelBulk_5_single_precision.h"

#include <domain_decomposition/BlockDataID.h>
#include <memory>

#include "utils.hpp"

namespace walberla {
namespace detail {
template <typename FloatType = double> struct KernelTrait {
  using ReactionKernelBulk_1 =
      pystencils::ReactionKernelBulk_1_double_precision;
  using ReactionKernelBulk_2 =
      pystencils::ReactionKernelBulk_2_double_precision;
  using ReactionKernelBulk_3 =
      pystencils::ReactionKernelBulk_3_double_precision;
  using ReactionKernelBulk_4 =
      pystencils::ReactionKernelBulk_4_double_precision;
  using ReactionKernelBulk_5 =
      pystencils::ReactionKernelBulk_5_double_precision;
};
template <> struct KernelTrait<float> {
  using ReactionKernelBulk_1 =
      pystencils::ReactionKernelBulk_1_single_precision;
  using ReactionKernelBulk_2 =
      pystencils::ReactionKernelBulk_2_single_precision;
  using ReactionKernelBulk_3 =
      pystencils::ReactionKernelBulk_3_single_precision;
  using ReactionKernelBulk_4 =
      pystencils::ReactionKernelBulk_4_single_precision;
  using ReactionKernelBulk_5 =
      pystencils::ReactionKernelBulk_5_single_precision;
};

template <typename FloatType>
auto get_kernel(const std::vector<std::shared_ptr<EKReactant>> &reactants,
                const double coefficient) {

  const auto coeff_casted = numeric_cast<FloatType>(coefficient);

  switch (reactants.size()) {
  case 1: {
    const auto &reactant = reactants[0];
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactant);

    auto kernel =
        std::make_shared<typename KernelTrait<FloatType>::ReactionKernelBulk_1>(
            density_id_0, order_0, coeff_casted, stoech_coeff_0);

    return kernel->getSweep();
  }
  case 2: {
    const auto [density_id_0, order_0, stoech_coeff_0] =
        detail::get_reaction_details<FloatType>(reactants[0]);
    const auto [density_id_1, order_1, stoech_coeff_1] =
        detail::get_reaction_details<FloatType>(reactants[1]);

    auto kernel =
        std::make_shared<typename KernelTrait<FloatType>::ReactionKernelBulk_2>(
            density_id_0, density_id_1, order_0, order_1, coeff_casted,
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

    auto kernel =
        std::make_shared<typename KernelTrait<FloatType>::ReactionKernelBulk_3>(
            density_id_0, density_id_1, density_id_2, order_0, order_1, order_2,
            coeff_casted, stoech_coeff_0, stoech_coeff_1, stoech_coeff_2);

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

    auto kernel =
        std::make_shared<typename KernelTrait<FloatType>::ReactionKernelBulk_4>(
            density_id_0, density_id_1, density_id_2, density_id_3, order_0,
            order_1, order_2, order_3, coeff_casted, stoech_coeff_0,
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

    auto kernel =
        std::make_shared<typename KernelTrait<FloatType>::ReactionKernelBulk_5>(
            density_id_0, density_id_1, density_id_2, density_id_3,
            density_id_4, order_0, order_1, order_2, order_3, order_4,
            coeff_casted, stoech_coeff_0, stoech_coeff_1, stoech_coeff_2,
            stoech_coeff_3, stoech_coeff_4);

    return kernel->getSweep();
  }
  default:
    throw std::runtime_error("reactions of this size are not implemented!");
  }
}

auto get_kernel_helper(
    const std::vector<std::shared_ptr<EKReactant>> &reactants,
    const double coefficient) {

  const auto is_double_precision =
      reactants[0]->get_species()->is_double_precision();

  if (is_double_precision) {
    return get_kernel<double>(reactants, coefficient);
  }

  return get_kernel<float>(reactants, coefficient);
}
} // namespace detail

void EKReactionImplBulk::perform_reaction() {
  // TODO: if my understanding is correct:
  //  the kernels need to either run in the ghost layers and do the
  //  synchronization before or not run and do a synchronization afterwards.
  //  The better solution is probably the latter one. Not sure why it fails
  //  atm.

  auto kernel = detail::get_kernel_helper(get_reactants(), get_coefficient());
  for (auto &block : *get_lattice()->get_blocks()) {
    kernel(&block);
  }
}
} // namespace walberla
