/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

/**
 * @file
 * Selection of the generated reaction kernels (see
 * @c maintainer/walberla_kernels/reactions.py) at runtime: the kernels exist
 * for 1 to 5 reactants, in single and double precision, for bulk reactions
 * (all cells of the interior of a block) and indexed reactions (the cells of
 * a @c walberla::v8::sweep::SparseIndexList), on the CPU and on the GPU.
 * The GPU classes are only available when the file is compiled by a CUDA
 * compiler.
 */

#include <waLBerlaDefinitions.h>

#include "../generated_kernels/ReactionKernelBulk1DoublePrecision.h"
#include "../generated_kernels/ReactionKernelBulk1SinglePrecision.h"
#include "../generated_kernels/ReactionKernelBulk2DoublePrecision.h"
#include "../generated_kernels/ReactionKernelBulk2SinglePrecision.h"
#include "../generated_kernels/ReactionKernelBulk3DoublePrecision.h"
#include "../generated_kernels/ReactionKernelBulk3SinglePrecision.h"
#include "../generated_kernels/ReactionKernelBulk4DoublePrecision.h"
#include "../generated_kernels/ReactionKernelBulk4SinglePrecision.h"
#include "../generated_kernels/ReactionKernelBulk5DoublePrecision.h"
#include "../generated_kernels/ReactionKernelBulk5SinglePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed1DoublePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed1SinglePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed2DoublePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed2SinglePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed3DoublePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed3SinglePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed4DoublePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed4SinglePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed5DoublePrecision.h"
#include "../generated_kernels/ReactionKernelIndexed5SinglePrecision.h"

#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
#include "../generated_kernels/ReactionKernelBulk1DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk1SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk2DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk2SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk3DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk3SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk4DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk4SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk5DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelBulk5SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed1DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed1SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed2DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed2SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed3DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed3SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed4DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed4SinglePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed5DoublePrecisionCUDA.h"
#include "../generated_kernels/ReactionKernelIndexed5SinglePrecisionCUDA.h"

#endif

#include <walberla_bridge/Architecture.hpp>

#include <core/DataTypes.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>

#include <walberla/v8/Memory.hpp>
#include <walberla/v8/sweep/SparseIndexList.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace walberla {
namespace detail {
namespace ReactionKernelSelector {

/** @brief The numbers of reactants that have a kernel. */
using Sizes = std::index_sequence<1, 2, 3, 4, 5>;
/** @brief Number of reactants of the largest reaction that has a kernel. */
constexpr std::size_t max_num_reactants = Sizes::size();

/**
 * @brief Index list of the indexed reaction kernels.
 * The index list of the GPU kernels is in unified memory, which the host
 * writes to and the kernel reads from without a copy.
 */
template <lbmpy::Arch Architecture>
using IndexList = walberla::v8::sweep::SparseIndexList<
    walberla::v8::sweep::CellIdx,
    std::conditional_t<Architecture == lbmpy::Arch::CPU,
                       walberla::v8::memory::memtag::host,
                       walberla::v8::memory::memtag::unified>>;

/** @brief Dependent @c false, to make a @c static_assert fail on use only. */
template <typename FloatType, lbmpy::Arch Architecture, std::size_t N>
inline constexpr bool always_false = false;

/**
 * @brief The kernels for @c N reactants.
 * Specialized for the floating-point type, the architecture and the number
 * of reactants. @c Bulk and @c Indexed are the kernels of bulk and of indexed
 * reactions. All kernels are constructed from the density fields of the
 * reactants (preceded by the index list for the indexed reactions), followed
 * by their orders, the rate coefficient and their stoichiometric
 * coefficients.
 */
template <typename FloatType, lbmpy::Arch Architecture, std::size_t N>
struct KernelTrait {
  static_assert(always_false<FloatType, Architecture, N>,
                "the GPU reaction kernels need a CUDA compiler (compile this "
                "translation unit with nvcc)");
};

template <> struct KernelTrait<double, lbmpy::Arch::CPU, 1> {
  using Bulk = pystencils::ReactionKernelBulk1DoublePrecision;
  using Indexed = pystencils::ReactionKernelIndexed1DoublePrecision;
};

template <> struct KernelTrait<double, lbmpy::Arch::CPU, 2> {
  using Bulk = pystencils::ReactionKernelBulk2DoublePrecision;
  using Indexed = pystencils::ReactionKernelIndexed2DoublePrecision;
};

template <> struct KernelTrait<double, lbmpy::Arch::CPU, 3> {
  using Bulk = pystencils::ReactionKernelBulk3DoublePrecision;
  using Indexed = pystencils::ReactionKernelIndexed3DoublePrecision;
};

template <> struct KernelTrait<double, lbmpy::Arch::CPU, 4> {
  using Bulk = pystencils::ReactionKernelBulk4DoublePrecision;
  using Indexed = pystencils::ReactionKernelIndexed4DoublePrecision;
};

template <> struct KernelTrait<double, lbmpy::Arch::CPU, 5> {
  using Bulk = pystencils::ReactionKernelBulk5DoublePrecision;
  using Indexed = pystencils::ReactionKernelIndexed5DoublePrecision;
};

template <> struct KernelTrait<float, lbmpy::Arch::CPU, 1> {
  using Bulk = pystencils::ReactionKernelBulk1SinglePrecision;
  using Indexed = pystencils::ReactionKernelIndexed1SinglePrecision;
};

template <> struct KernelTrait<float, lbmpy::Arch::CPU, 2> {
  using Bulk = pystencils::ReactionKernelBulk2SinglePrecision;
  using Indexed = pystencils::ReactionKernelIndexed2SinglePrecision;
};

template <> struct KernelTrait<float, lbmpy::Arch::CPU, 3> {
  using Bulk = pystencils::ReactionKernelBulk3SinglePrecision;
  using Indexed = pystencils::ReactionKernelIndexed3SinglePrecision;
};

template <> struct KernelTrait<float, lbmpy::Arch::CPU, 4> {
  using Bulk = pystencils::ReactionKernelBulk4SinglePrecision;
  using Indexed = pystencils::ReactionKernelIndexed4SinglePrecision;
};

template <> struct KernelTrait<float, lbmpy::Arch::CPU, 5> {
  using Bulk = pystencils::ReactionKernelBulk5SinglePrecision;
  using Indexed = pystencils::ReactionKernelIndexed5SinglePrecision;
};

#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)

template <> struct KernelTrait<double, lbmpy::Arch::GPU, 1> {
  using Bulk = pystencils::ReactionKernelBulk1DoublePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed1DoublePrecisionCUDA;
};

template <> struct KernelTrait<double, lbmpy::Arch::GPU, 2> {
  using Bulk = pystencils::ReactionKernelBulk2DoublePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed2DoublePrecisionCUDA;
};

template <> struct KernelTrait<double, lbmpy::Arch::GPU, 3> {
  using Bulk = pystencils::ReactionKernelBulk3DoublePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed3DoublePrecisionCUDA;
};

template <> struct KernelTrait<double, lbmpy::Arch::GPU, 4> {
  using Bulk = pystencils::ReactionKernelBulk4DoublePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed4DoublePrecisionCUDA;
};

template <> struct KernelTrait<double, lbmpy::Arch::GPU, 5> {
  using Bulk = pystencils::ReactionKernelBulk5DoublePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed5DoublePrecisionCUDA;
};

template <> struct KernelTrait<float, lbmpy::Arch::GPU, 1> {
  using Bulk = pystencils::ReactionKernelBulk1SinglePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed1SinglePrecisionCUDA;
};

template <> struct KernelTrait<float, lbmpy::Arch::GPU, 2> {
  using Bulk = pystencils::ReactionKernelBulk2SinglePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed2SinglePrecisionCUDA;
};

template <> struct KernelTrait<float, lbmpy::Arch::GPU, 3> {
  using Bulk = pystencils::ReactionKernelBulk3SinglePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed3SinglePrecisionCUDA;
};

template <> struct KernelTrait<float, lbmpy::Arch::GPU, 4> {
  using Bulk = pystencils::ReactionKernelBulk4SinglePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed4SinglePrecisionCUDA;
};

template <> struct KernelTrait<float, lbmpy::Arch::GPU, 5> {
  using Bulk = pystencils::ReactionKernelBulk5SinglePrecisionCUDA;
  using Indexed = pystencils::ReactionKernelIndexed5SinglePrecisionCUDA;
};

#endif

/**
 * @brief Create the kernel of the type @c Kernel for the reactants.
 * @param reactants Reactants of the reaction, @c Is... their indices
 * @param coefficient Rate coefficient of the reaction
 * @param leading_args Leading arguments of the kernel (the index list)
 * @return the sweep that runs the kernel on a block
 */
template <typename Kernel, typename FloatType, class Reactant,
          std::size_t... Is, typename... LeadingArgs>
std::function<void(IBlock *)>
make_sweep(std::vector<std::shared_ptr<Reactant>> const &reactants,
           double const coefficient, std::index_sequence<Is...>,
           LeadingArgs const &...leading_args) {
  auto kernel = std::make_shared<Kernel>(
      leading_args...,
      walberla::BlockDataID(reactants[Is]->get_species()->get_density_id())...,
      numeric_cast<FloatType>(reactants[Is]->get_order())...,
      numeric_cast<FloatType>(coefficient),
      numeric_cast<FloatType>(reactants[Is]->get_stoech_coeff())...);
  return [kernel](IBlock *block) { (*kernel)(block); };
}

/**
 * @brief Select the kernel for the number of reactants.
 * @param size Number of reactants, which is one of @c Ns...
 * @param build Function object whose
 *        <tt>operator()<FloatType, N>()</tt> creates the sweep for @c N
 *        reactants
 */
template <typename FloatType, class Build, std::size_t... Ns>
std::function<void(IBlock *)> select_size(std::size_t const size, Build &build,
                                          std::index_sequence<Ns...>) {
  std::function<void(IBlock *)> sweep;
  ((size == Ns ? void(sweep = build.template operator()<FloatType, Ns>())
               : void()),
   ...);
  return sweep;
}

/**
 * @brief Select the kernel for the floating-point type and the number of
 * reactants of a reaction.
 * @throw std::runtime_error if the reaction has no reactants or more than
 *        @ref max_num_reactants
 */
template <class Reactant, class Build>
std::function<void(IBlock *)>
select_kernel(std::vector<std::shared_ptr<Reactant>> const &reactants,
              Build &&build) {
  auto const size = reactants.size();
  if (size == 0u or size > max_num_reactants) {
    throw std::runtime_error("reactions of this size are not implemented!");
  }
  if (reactants[0]->get_species()->is_double_precision()) {
    return select_size<double>(size, build, Sizes{});
  }
  return select_size<float>(size, build, Sizes{});
}

/**
 * @brief Get the sweep of a bulk reaction, which runs on all cells of the
 * interior of a block.
 * The floating-point type is the one of the first reactant.
 */
template <lbmpy::Arch Architecture, class Reactant>
std::function<void(IBlock *)>
get_bulk_kernel(std::vector<std::shared_ptr<Reactant>> const &reactants,
                double const coefficient) {
  return select_kernel(reactants, [&]<typename FloatType, std::size_t N>() {
    using Kernel = typename KernelTrait<FloatType, Architecture, N>::Bulk;
    return make_sweep<Kernel, FloatType>(reactants, coefficient,
                                         std::make_index_sequence<N>{});
  });
}

/**
 * @brief Get the sweep of an indexed reaction, which runs on the cells of
 * an index list.
 * The floating-point type is the one of the first reactant.
 * The index list is copied by handle: a @c SparseIndexList only holds the
 * @c BlockDataID of block data that is owned by the block forest.
 */
template <lbmpy::Arch Architecture, class Reactant>
std::function<void(IBlock *)>
get_indexed_kernel(std::vector<std::shared_ptr<Reactant>> const &reactants,
                   double const coefficient,
                   IndexList<Architecture> const &index_list) {
  return select_kernel(reactants, [&]<typename FloatType, std::size_t N>() {
    using Kernel = typename KernelTrait<FloatType, Architecture, N>::Indexed;
    return make_sweep<Kernel, FloatType>(
        reactants, coefficient, std::make_index_sequence<N>{}, index_list);
  });
}

} // namespace ReactionKernelSelector
} // namespace detail
} // namespace walberla
