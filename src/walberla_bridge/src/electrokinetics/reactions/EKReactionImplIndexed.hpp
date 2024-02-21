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

#include "generated_kernels/ReactionKernelIndexed_all.h"

#include <walberla_bridge/BlockAndCell.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactant.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBaseIndexed.hpp>

#include <blockforest/StructuredBlockForest.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/AddToStorage.h>
#include <field/FlagField.h>
#include <field/FlagUID.h>

#include <utils/Vector.hpp>

#include <cassert>
#include <cstddef>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

namespace walberla {

class EKReactionImplIndexed : public EKReactionBaseIndexed {
private:
  BlockDataID m_flagfield_id;
  BlockDataID m_indexvector_id;
  bool m_pending_changes;

public:
  /** Flag for domain cells, i.e. all cells. */
  FlagUID const Domain_flag{"domain"};
  /** Flag for boundary cells. */
  FlagUID const Boundary_flag{"boundary"};

  using FlagField = field::FlagField<uint8_t>;
  using IndexVectors = detail::ReactionKernelIndexedSelector::KernelTrait<>::
      ReactionKernelIndexed::IndexVectors;
  using IndexInfo = detail::ReactionKernelIndexedSelector::KernelTrait<>::
      ReactionKernelIndexed::IndexInfo;

private:
  auto get_flag_field_and_flag(IBlock *block, BlockDataID const &flagfield_id) {
    auto const flag_field =
        block->template uncheckedFastGetData<FlagField>(flagfield_id);
    auto const boundary_flag = flag_field->getFlag(Boundary_flag);
    return std::make_tuple(flag_field, boundary_flag);
  }

public:
  EKReactionImplIndexed(std::shared_ptr<LatticeWalberla> const &lattice,
                        reactants_type const &reactants, double coefficient)
      : EKReactionBaseIndexed(lattice, reactants, coefficient),
        m_pending_changes(false) {
    m_flagfield_id = field::addFlagFieldToStorage<FlagField>(
        get_lattice()->get_blocks(), "flag field reaction",
        get_lattice()->get_ghost_layers());

    auto createIdxVector = [](IBlock *const, StructuredBlockStorage *const) {
      return new IndexVectors();
    };
    m_indexvector_id = get_lattice()
                           ->get_blocks()
                           ->template addStructuredBlockData<IndexVectors>(
                               createIdxVector, "IndexField");

    for (auto &block : *get_lattice()->get_blocks()) {
      auto flag_field = block.template getData<FlagField>(m_flagfield_id);
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
  ~EKReactionImplIndexed() override = default;

  using EKReactionBaseIndexed::get_coefficient;
  using EKReactionBaseIndexed::get_lattice;
  using EKReactionBaseIndexed::get_reactants;

  void perform_reaction() override {
    boundary_update();
    auto kernel = detail::ReactionKernelIndexedSelector::get_kernel(
        get_reactants(), get_coefficient(), m_indexvector_id);
    for (auto &block : *get_lattice()->get_blocks()) {
      kernel(&block);
    }
  }

  void set_node_is_boundary(Utils::Vector3i const &node,
                            bool is_boundary) override {
    if (auto bc = get_block_and_cell(*get_lattice(), node, true)) {
      auto const [flag_field, boundary_flag] =
          get_flag_field_and_flag(bc->block, m_flagfield_id);
      if (is_boundary) {
        flag_field->addFlag(bc->cell, boundary_flag);
      } else {
        flag_field->removeFlag(bc->cell, boundary_flag);
      }
      m_pending_changes = true;
    }
  }

  [[nodiscard]] std::optional<bool>
  get_node_is_boundary(Utils::Vector3i const &node) override {
    if (auto bc = get_block_and_cell(*get_lattice(), node, true)) {
      auto const [flag_field, boundary_flag] =
          get_flag_field_and_flag(bc->block, m_flagfield_id);
      return {flag_field->isFlagSet(bc->cell, boundary_flag)};
    }
    return std::nullopt;
  }

  void boundary_update() {
    if (m_pending_changes) {
      for (auto &block : *get_lattice()->get_blocks()) {
        fillFromFlagField(block);
      }
      m_pending_changes = false;
    }
  }

private:
  void fillFromFlagField(IBlock &block) {
    auto *indexVectors =
        block.uncheckedFastGetData<IndexVectors>(m_indexvector_id);
    auto &indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
    auto &indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
    auto &indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);

    auto *flagField = block.getData<FlagField>(m_flagfield_id);

    assert(flagField->flagExists(Boundary_flag) and
           flagField->flagExists(Domain_flag));

    auto boundaryFlag = flagField->getFlag(Boundary_flag);
    auto domainFlag = flagField->getFlag(Domain_flag);

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
        if (inner.contains(it.x(), it.y(), it.z())) {
          indexVectorInner.push_back(element);
        } else {
          indexVectorOuter.push_back(element);
        }
      }
    }

    indexVectors->syncGPU();
  }
};

} // namespace walberla
