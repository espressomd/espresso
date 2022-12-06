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

#include "EKReactionImplIndexed.hpp"

#include "walberla_bridge/BlockAndCell.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"
#include "walberla_bridge/electrokinetics/reactions/EKReactant.hpp"
#include "walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp"

#include "generated_kernels/ReactionKernelIndexed_all.h"

#include <domain_decomposition/BlockDataID.h>
#include <field/AddToStorage.h>

#include <boost/optional.hpp>

#include <memory>
#include <tuple>
#include <vector>

namespace walberla {

/// Flag for domain cells, i.e. all cells
const FlagUID Domain_flag("domain");
/// Flag for boundary cells
const FlagUID Boundary_flag("boundary");

namespace detail {
// FlagField to use
using FlagField = FlagField<uint8_t>;

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

EKReactionImplIndexed::EKReactionImplIndexed(
    std::shared_ptr<LatticeWalberla> lattice,
    std::vector<std::shared_ptr<EKReactant>> reactants, double coefficient)
    : EKReactionBase(lattice, reactants, coefficient),
      m_pending_changes(false) {
  m_flagfield_id =
      static_cast<std::size_t>(field::addFlagFieldToStorage<detail::FlagField>(
          get_lattice()->get_blocks(), "flag field reaction",
          get_lattice()->get_ghost_layers()));

  // take one IndexVector as a dummy-value
  using IndexVectors = detail::ReactionKernelIndexedSelector::KernelTrait<>::
      ReactionKernelIndexed::IndexVectors;

  auto createIdxVector = [](IBlock *const, StructuredBlockStorage *const) {
    return new IndexVectors();
  };
  m_indexvector_id = static_cast<std::size_t>(
      get_lattice()
          ->get_blocks()
          ->template addStructuredBlockData<IndexVectors>(createIdxVector,
                                                          "IndexField"));

  for (auto &block : *get_lattice()->get_blocks()) {
    auto flag_field =
        block.template getData<detail::FlagField>(BlockDataID(m_flagfield_id));
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

void EKReactionImplIndexed::perform_reaction() {
  boundary_update();

  auto kernel = detail::ReactionKernelIndexedSelector::get_kernel(
      get_reactants(), get_coefficient(), BlockDataID(get_indexvector_id()));

  for (auto &block : *get_lattice()->get_blocks()) {
    kernel(&block);
  }
}

void EKReactionImplIndexed::set_node_is_boundary(const Utils::Vector3i &node,
                                                 bool is_boundary) {
  auto bc = get_block_and_cell(*get_lattice(), node, true);
  if (!bc)
    return;

  auto [flag_field, boundary_flag] =
      detail::get_flag_field_and_flag<detail::FlagField>(
          bc->block, BlockDataID(get_flagfield_id()));
  if (is_boundary) {
    flag_field->addFlag(bc->cell, boundary_flag);
  } else {
    flag_field->removeFlag(bc->cell, boundary_flag);
  }
  m_pending_changes = true;
}

boost::optional<bool>
EKReactionImplIndexed::get_node_is_boundary(const Utils::Vector3i &node) {
  auto bc = get_block_and_cell(*get_lattice(), node, true);
  if (!bc)
    return {boost::none};

  auto [flag_field, boundary_flag] =
      detail::get_flag_field_and_flag<detail::FlagField>(
          bc->block, BlockDataID(get_flagfield_id()));
  return {flag_field->isFlagSet(bc->cell, boundary_flag)};
}

void EKReactionImplIndexed::boundary_update() {
  // take one IndexVector/IndexInfo as a dummy-value
  using IndexVectors = detail::ReactionKernelIndexedSelector::KernelTrait<>::
      ReactionKernelIndexed::IndexVectors;
  using IndexInfo = detail::ReactionKernelIndexedSelector::KernelTrait<>::
      ReactionKernelIndexed::IndexInfo;

  if (m_pending_changes) {
    detail::fillFromFlagField<detail::FlagField, IndexVectors, IndexInfo>(
        get_lattice()->get_blocks(), BlockDataID(get_indexvector_id()),
        BlockDataID(get_flagfield_id()), Boundary_flag, Domain_flag);
    m_pending_changes = false;
  }
}
} // namespace walberla
