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

// kernel generated with pystencils v1.2, lbmpy v1.2,
// lbmpy_walberla/pystencils_walberla from waLBerla commit ref:
// a839fac6ef7d0c58e7710e4d50490e9dd7146b4a

#pragma once

#include "ReactionKernelBulk_1_double_precision.h"
#include "ReactionKernelBulk_1_single_precision.h"

#include "ReactionKernelBulk_2_double_precision.h"
#include "ReactionKernelBulk_2_single_precision.h"

#include "ReactionKernelBulk_3_double_precision.h"
#include "ReactionKernelBulk_3_single_precision.h"

#include "ReactionKernelBulk_4_double_precision.h"
#include "ReactionKernelBulk_4_single_precision.h"

#include "ReactionKernelBulk_5_double_precision.h"
#include "ReactionKernelBulk_5_single_precision.h"

#include <domain_decomposition/BlockDataID.h>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace walberla {
namespace detail {
namespace ReactionKernelBulkSelector {

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

template <typename FloatType, class Reactant, std::size_t... ints>
auto get_kernel_impl(const std::vector<std::shared_ptr<Reactant>> &reactants,
                     const double coefficient,
                     std::index_sequence<ints...> int_seq) {
  auto kernel = std::make_shared<
      typename KernelTrait<FloatType, int_seq.size()>::ReactionKernelBulk>(
      walberla::BlockDataID(
          reactants[ints]->get_species()->get_density_id())...,
      numeric_cast<FloatType>(reactants[ints]->get_order())...,
      numeric_cast<FloatType>(coefficient),
      numeric_cast<FloatType>(reactants[ints]->get_stoech_coeff())...);

  std::function<void(IBlock *)> sweep = [kernel](IBlock *b) { kernel->run(b); };
  return sweep;
}

template <typename FloatType, class Reactant, class... Args>
auto get_kernel_impl(const std::vector<std::shared_ptr<Reactant>> &reactants,
                     Args... args) {
  switch (reactants.size()) {

  case 1:
    return get_kernel_impl<FloatType>(reactants, args...,
                                      std::make_index_sequence<1>{});

  case 2:
    return get_kernel_impl<FloatType>(reactants, args...,
                                      std::make_index_sequence<2>{});

  case 3:
    return get_kernel_impl<FloatType>(reactants, args...,
                                      std::make_index_sequence<3>{});

  case 4:
    return get_kernel_impl<FloatType>(reactants, args...,
                                      std::make_index_sequence<4>{});

  case 5:
    return get_kernel_impl<FloatType>(reactants, args...,
                                      std::make_index_sequence<5>{});

  default:
    throw std::runtime_error("reactions of this size are not implemented!");
  }
}

template <class Reactant, class... Args>
auto get_kernel(const std::vector<std::shared_ptr<Reactant>> &reactants,
                Args... args) {

  const auto is_double_precision =
      reactants[0]->get_species()->is_double_precision();

  if (is_double_precision) {
    return get_kernel_impl<double>(reactants, args...);
  }

  return get_kernel_impl<float>(reactants, args...);
}

} // namespace ReactionKernelBulkSelector
} // namespace detail
} // namespace walberla