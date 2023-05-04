/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include <walberla_bridge/BlockAndCell.hpp>
#include <walberla_bridge/utils/walberla_utils.hpp>

#include <blockforest/StructuredBlockForest.h>
#include <field/FlagField.h>

#include <utils/Vector.hpp>

#include <cassert>
#include <functional>
#include <memory>
#include <tuple>
#include <unordered_map>

namespace walberla {

/// Flag for domain cells, i.e. all cells
FlagUID const Domain_flag("domain");
/// Flag for boundary cells
FlagUID const Boundary_flag("boundary");

template <typename T, typename BoundaryClass> class BoundaryHandling {
  /** Container for the map between cells and values. */
  class DynamicValueCallback {
  public:
    DynamicValueCallback() {
      m_value_boundary = std::make_shared<std::unordered_map<Cell, T>>();
    }

    [[nodiscard]] T operator()(
        Cell const &local,
        std::shared_ptr<blockforest::StructuredBlockForest> const &blocks,
        IBlock &block) const {
      Cell global;
      blocks->transformBlockLocalToGlobalCell(global, block, local);
      return get_value(global);
    }

    template <typename U>
    void set_node_boundary_value(Utils::Vector3i const &node, U const &val) {
      auto const global = Cell(node[0], node[1], node[2]);
      (*m_value_boundary)[global] = es2walberla<U, T>(val);
    }

    void unset_node_boundary_value(Utils::Vector3i const &node) {
      auto const global = Cell(node[0], node[1], node[2]);
      assert(m_value_boundary->count(global));
      m_value_boundary->erase(global);
    }

    [[nodiscard]] auto
    get_node_boundary_value(Utils::Vector3i const &node) const {
      auto const global = Cell(node[0], node[1], node[2]);
      return walberla2es(get_value(global));
    }

    bool node_is_boundary(Utils::Vector3i const &node) const {
      auto const global = Cell(node[0], node[1], node[2]);
      return m_value_boundary->count(global) != 0;
    }

  private:
    std::shared_ptr<std::unordered_map<Cell, T>> m_value_boundary;
    static constexpr T default_value{};

    [[nodiscard]] T get_value(Cell const &cell) const {
      if (m_value_boundary->count(cell) == 0) {
        return default_value;
      }
      return m_value_boundary->at(cell);
    }
  };

  [[nodiscard]] inline auto get_flag_field_and_flag(IBlock *block) const {
    auto const flag_field =
        block->template uncheckedFastGetData<FlagField>(m_flag_field_id);
    auto const boundary_flag = flag_field->getFlag(Boundary_flag);
    return std::make_tuple(flag_field, boundary_flag);
  }

public:
  using value_type = T;
  using FlagField = field::FlagField<uint8_t>;

  BoundaryHandling(std::shared_ptr<StructuredBlockForest> blocks,
                   BlockDataID value_field_id, BlockDataID flag_field_id)
      : m_blocks(std::move(blocks)), m_value_field_id(value_field_id),
        m_flag_field_id(flag_field_id), m_callback(DynamicValueCallback()),
        m_pending_changes(false) {
    // reinitialize the flag field
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b) {
      flag_reset_kernel(&*b);
    }
    // instantiate the boundary sweep
    std::function callback = m_callback;
    m_boundary =
        std::make_shared<BoundaryClass>(m_blocks, m_value_field_id, callback);
  }

  void operator()(IBlock *block) { (*m_boundary)(block); }

  [[nodiscard]] bool node_is_boundary(Utils::Vector3i const &node) const {
    return m_callback.node_is_boundary(node);
  }

  [[nodiscard]] auto
  get_node_value_at_boundary(Utils::Vector3i const &node) const {
    return m_callback.get_node_boundary_value(node);
  }

  template <typename U>
  void set_node_value_at_boundary(Utils::Vector3i const &node, U const &v,
                                  BlockAndCell const &bc) {
    auto [flag_field, boundary_flag] = get_flag_field_and_flag(bc.block);
    m_callback.set_node_boundary_value(node, v);
    flag_field->addFlag(bc.cell, boundary_flag);
    m_pending_changes = true;
  }

  void remove_node_from_boundary(Utils::Vector3i const &node,
                                 BlockAndCell const &bc) {
    auto [flag_field, boundary_flag] = get_flag_field_and_flag(bc.block);
    m_callback.unset_node_boundary_value(node);
    flag_field->removeFlag(bc.cell, boundary_flag);
    m_pending_changes = true;
  }

  /** Assign boundary conditions to boundary cells. */
  void boundary_update() {
    if (m_pending_changes) {
      m_boundary->template fillFromFlagField<FlagField>(
          m_blocks, m_flag_field_id, Boundary_flag, Domain_flag);
      m_pending_changes = false;
    }
  }

  std::tuple<int, int, int> block_dims(IBlock const &block) const {
    auto const field = block.template getData<FlagField>(m_flag_field_id);
    return {static_cast<int>(field->xSize()), static_cast<int>(field->ySize()),
            static_cast<int>(field->zSize())};
  }

private:
  std::shared_ptr<StructuredBlockForest> m_blocks;
  BlockDataID m_value_field_id;
  BlockDataID m_flag_field_id;
  DynamicValueCallback m_callback;
  std::shared_ptr<BoundaryClass> m_boundary;
  bool m_pending_changes;

  /** Register flags and set all cells to @ref Domain_flag. */
  void flag_reset_kernel(IBlock *const block) {
    auto flag_field = block->template getData<FlagField>(m_flag_field_id);
    // register flags
    if (!flag_field->flagExists(Domain_flag))
      flag_field->registerFlag(Domain_flag);
    if (!flag_field->flagExists(Boundary_flag))
      flag_field->registerFlag(Boundary_flag);
    // mark all cells as domain cells and fluid cells
    auto domain_flag = flag_field->getFlag(Domain_flag);
    auto boundary_flag = flag_field->getFlag(Boundary_flag);
    for (auto it = flag_field->begin(); it != flag_field->end(); ++it) {
      flag_field->addFlag(it.x(), it.y(), it.z(), domain_flag);
      flag_field->removeFlag(it.x(), it.y(), it.z(), boundary_flag);
    }
  }
};

} // namespace walberla
