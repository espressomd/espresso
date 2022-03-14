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
#include "BlockAndCell.hpp"
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

#include <domain_decomposition/BlockDataID.h>
#include <field/AddToStorage.h>

namespace walberla {

/// Flag for domain cells, i.e. all cells
const FlagUID Domain_flag("domain");
/// Flag for boundary cells
const FlagUID Boundary_flag("boundary");

namespace detail {
// FlagField to use
using FlagField = FlagField<uint8_t>;

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

template <typename FlagField>
inline auto
get_flag_field_and_flag(IBlock *block,
                        const domain_decomposition::BlockDataID &flagfield_id) {
  auto const flag_field =
      block->template uncheckedFastGetData<FlagField>(flagfield_id);
  auto const boundary_flag = flag_field->getFlag(Boundary_flag);
  return std::make_tuple(flag_field, boundary_flag);
}

template <typename FlagField, typename IndexVectors, typename IndexInfo>
void fillFromFlagField(IBlock *block, BlockDataID indexVectorID,
                       ConstBlockDataID flagFieldID, FlagUID boundaryFlagUID,
                       FlagUID domainFlagUID) {
  auto *indexVectors = block->uncheckedFastGetData<IndexVectors>(indexVectorID);
  auto &indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
  auto &indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
  auto &indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);

  auto *flagField = block->getData<FlagField>(flagFieldID);

  if (!(flagField->flagExists(boundaryFlagUID) &&
        flagField->flagExists(domainFlagUID)))
    return;

  auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
  auto domainFlag = flagField->getFlag(domainFlagUID);

  auto inner = flagField->xyzSize();
  inner.expand(cell_idx_t(-1));

  indexVectorAll.clear();
  indexVectorInner.clear();
  indexVectorOuter.clear();

  auto flagWithGLayers = flagField->xyzSizeWithGhostLayer();
  for (auto it = flagField->beginWithGhostLayerXYZ(); it != flagField->end();
       ++it) {

    if (!isFlagSet(it, boundaryFlag))
      continue;
    if (flagWithGLayers.contains(it.x(), it.y(), it.z()) &&
        isFlagSet(it.neighbor(0, 0, 0, 0), domainFlag)) {

      auto element = IndexInfo(it.x(), it.y(), it.z());

      indexVectorAll.push_back(element);
      if (inner.contains(it.x(), it.y(), it.z()))
        indexVectorInner.push_back(element);
      else
        indexVectorOuter.push_back(element);
    }
  }

  indexVectors->syncGPU();
}

template <typename FlagField, typename IndexVectors, typename IndexInfo>
void fillFromFlagField(const shared_ptr<StructuredBlockForest> &blocks,
                       BlockDataID indexVectorID, ConstBlockDataID flagFieldID,
                       FlagUID boundaryFlagUID, FlagUID domainFlagUID) {
  for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
    fillFromFlagField<FlagField, IndexVectors, IndexInfo>(
        blockIt.get(), indexVectorID, flagFieldID, boundaryFlagUID,
        domainFlagUID);
}
} // namespace detail

template <typename FloatType>
EKReactionIndexed<FloatType>::EKReactionIndexed(
    std::shared_ptr<LatticeWalberla> lattice,
    std::vector<std::shared_ptr<EKReactant<FloatType>>> reactants,
    FloatType coefficient)
    : EKReactionBase<FloatType>(lattice, reactants, coefficient),
      m_pending_changes(false) {
  m_flagfield_id = field::addFlagFieldToStorage<detail::FlagField>(
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

  for (auto &block : *get_lattice()->get_blocks()) {
    auto flag_field = block.template getData<detail::FlagField>(m_flagfield_id);
    // register flags
    flag_field->registerFlag(Domain_flag);
    flag_field->registerFlag(Boundary_flag);
    // mark all cells as domain cells and fluid cells
    auto domain_flag = flag_field->getFlag(Domain_flag);
    auto boundary_flag = flag_field->getFlag(Boundary_flag);
    for (auto it = flag_field->begin(); it != flag_field->end(); ++it) {
      flag_field->addFlag(it.x(), it.y(), it.z(), domain_flag);
      flag_field->removeFlag(it.x(), it.y(), it.z(), boundary_flag);
    }
  }
}

template <typename FloatType>
void EKReactionIndexed<FloatType>::perform_reaction() {
  boundary_update();

  auto kernel =
      detail::get_kernel(get_reactants(), get_coefficient(), m_indexvector_id);

  for (auto &block : *get_lattice()->get_blocks()) {
    kernel(&block);
  }
}

template <typename FloatType>
void EKReactionIndexed<FloatType>::set_node_is_boundary(
    const Utils::Vector3i &node, bool is_boundary) {
  auto bc = get_block_and_cell(*get_lattice(), node, true);
  if (!bc)
    return;

  auto [flag_field, boundary_flag] =
      detail::get_flag_field_and_flag<detail::FlagField>(bc->block,
                                                         get_flagfield_id());
  if (is_boundary) {
    flag_field->addFlag(bc->cell, boundary_flag);
  } else {
    flag_field->removeFlag(bc->cell, boundary_flag);
  }
  m_pending_changes = true;
}

template <typename FloatType>
boost::optional<bool> EKReactionIndexed<FloatType>::get_node_is_boundary(
    const Utils::Vector3i &node) {
  auto bc = get_block_and_cell(*get_lattice(), node, true);
  if (!bc)
    return {boost::none};

  auto [flag_field, boundary_flag] =
      detail::get_flag_field_and_flag<detail::FlagField>(bc->block,
                                                         get_flagfield_id());
  return {flag_field->isFlagSet(bc->cell, boundary_flag)};
}

template <typename FloatType>
void EKReactionIndexed<FloatType>::boundary_update() {
  using IndexVectors = pystencils::ReactionKernelIndexed_1::IndexVectors;
  using IndexInfo = pystencils::ReactionKernelIndexed_1::IndexInfo;

  if (m_pending_changes) {
    detail::fillFromFlagField<detail::FlagField, IndexVectors, IndexInfo>(
        get_lattice()->get_blocks(), get_indexvector_id(), get_flagfield_id(),
        Boundary_flag, Domain_flag);
    m_pending_changes = false;
  }
}
} // namespace walberla