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

#include <blockforest/StructuredBlockForest.h>

#include <cstddef>
#include <memory>
#include <utility>

namespace walberla {
namespace detail {
template <typename FloatType = double, std::size_t N = 1> struct KernelTrait {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_1_double_precision;
};
template <> struct KernelTrait<double, 2> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_2_double_precision;
};
template <> struct KernelTrait<double, 3> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_3_double_precision;
};
template <> struct KernelTrait<double, 4> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_4_double_precision;
};
template <> struct KernelTrait<double, 5> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_5_double_precision;
};
template <> struct KernelTrait<float, 1> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_1_single_precision;
};
template <> struct KernelTrait<float, 2> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_2_single_precision;
};
template <> struct KernelTrait<float, 3> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_3_single_precision;
};
template <> struct KernelTrait<float, 4> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_4_single_precision;
};
template <> struct KernelTrait<float, 5> {
  using ReactionKernelBulk = pystencils::ReactionKernelBulk_5_single_precision;
};

template <typename FloatType, std::size_t... ints>
auto get_kernel(std::vector<std::shared_ptr<EKReactant>> const &reactants,
                double coefficient, std::index_sequence<ints...> int_seq) {
  auto kernel = std::make_shared<
      typename KernelTrait<FloatType, int_seq.size()>::ReactionKernelBulk>(
      walberla::BlockDataID(
          reactants[ints]->get_species()->get_density_id())...,
      numeric_cast<FloatType>(reactants[ints]->get_order())...,
      numeric_cast<FloatType>(coefficient),
      numeric_cast<FloatType>(reactants[ints]->get_stoech_coeff())...);

  return kernel->getSweep();
}

template <typename FloatType>
auto get_kernel(const std::vector<std::shared_ptr<EKReactant>> &reactants,
                const double coefficient) {
  switch (reactants.size()) {
  case 1:
    return get_kernel<FloatType>(reactants, coefficient,
                                 std::make_index_sequence<1>{});
  case 2:
    return get_kernel<FloatType>(reactants, coefficient,
                                 std::make_index_sequence<2>{});
  case 3:
    return get_kernel<FloatType>(reactants, coefficient,
                                 std::make_index_sequence<3>{});
  case 4:
    return get_kernel<FloatType>(reactants, coefficient,
                                 std::make_index_sequence<4>{});
  case 5:
    return get_kernel<FloatType>(reactants, coefficient,
                                 std::make_index_sequence<5>{});
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
