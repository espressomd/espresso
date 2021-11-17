/*
 * Copyright (C) 2021 The ESPResSo project
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

#include <blockforest/StructuredBlockForest.h>
#include <field/FlagField.h>

#include "BlockAndCell.hpp"
#include "walberla_utils.hpp"

#include "generated_kernels/Dynamic_UBB_double_precision.h"
#include "generated_kernels/Dynamic_UBB_single_precision.h"

#include <utils/Vector.hpp>

#include <cassert>
#include <functional>
#include <memory>
#include <tuple>
#include <unordered_map>

namespace walberla {

/// Flag for domain cells, i.e. all cells
const FlagUID Domain_flag("domain");
/// Flag for boundary cells
const FlagUID Boundary_flag("velocity bounce back");

namespace detail {
template <typename FT> struct BoundaryHandlingTrait {
  using Dynamic_UBB = lbm::Dynamic_UBB_double_precision;
};
template <> struct BoundaryHandlingTrait<float> {
  using Dynamic_UBB = lbm::Dynamic_UBB_single_precision;
};
} // namespace detail

/** Velocity boundary conditions sweep. */
template <typename FloatType> class BoundaryHandling {

  using Dynamic_UBB =
      typename detail::BoundaryHandlingTrait<FloatType>::Dynamic_UBB;

  /** Container for the map between cells and velocities. */
  class DynamicVelocityCallback {
  public:
    DynamicVelocityCallback() {
      m_velocity_boundary =
          std::make_shared<std::unordered_map<Cell, Vector3<FloatType>>>();
    }

    Vector3<FloatType> operator()(
        Cell const &local,
        std::shared_ptr<blockforest::StructuredBlockForest> const &blocks,
        IBlock &block) const {
      Cell global;
      blocks->transformBlockLocalToGlobalCell(global, block, local);
      return get_velocity(global);
    }

    void set_node_boundary_velocity(Utils::Vector3i const &node,
                                    Utils::Vector3d const &vel) {
      auto const global = Cell(node[0], node[1], node[2]);
      (*m_velocity_boundary)[global] = to_vector3<FloatType>(vel);
    }

    void unset_node_boundary_velocity(Utils::Vector3i const &node) {
      auto const global = Cell(node[0], node[1], node[2]);
      assert(m_velocity_boundary->count(global));
      m_velocity_boundary->erase(global);
    }

    Utils::Vector3d
    get_node_boundary_velocity(Utils::Vector3i const &node) const {
      auto const global = Cell(node[0], node[1], node[2]);
      return to_vector3d(get_velocity(global));
    }

  private:
    std::shared_ptr<std::unordered_map<Cell, Vector3<FloatType>>>
        m_velocity_boundary;
    static constexpr Vector3<FloatType> no_slip{0, 0, 0};

    Vector3<FloatType> get_velocity(Cell const &cell) const {
      if (m_velocity_boundary->count(cell) == 0) {
        return no_slip;
      }
      return m_velocity_boundary->at(cell);
    }
  };

  inline auto get_flag_field_and_flag(IBlock *block) const {
    auto const flag_field = block->template getData<FlagField>(m_flag_field_id);
    auto const boundary_flag = flag_field->getFlag(Boundary_flag);
    return std::make_tuple(flag_field, boundary_flag);
  }

public:
  using FlagField = field::FlagField<uint8_t>;

  BoundaryHandling(std::shared_ptr<StructuredBlockForest> blocks,
                   BlockDataID pdf_field_id, BlockDataID flag_field_id)
      : m_blocks(std::move(blocks)), m_pdf_field_id(pdf_field_id),
        m_flag_field_id(flag_field_id), m_callback(DynamicVelocityCallback()) {
    // reinitialize the flag field
    for (auto b = m_blocks->begin(); b != m_blocks->end(); ++b) {
      flag_reset_kernel(&*b);
    }
    // instantiate the boundary sweep
    std::function callback = m_callback;
    m_boundary =
        std::make_shared<Dynamic_UBB>(m_blocks, m_pdf_field_id, callback);
  }

  void operator()(IBlock *block) { (*m_boundary)(block); }

  bool node_is_boundary(BlockAndCell const &bc) const {
    auto [flag_field, boundary_flag] = get_flag_field_and_flag(bc.block);
    return flag_field->isFlagSet(bc.cell, boundary_flag);
  }

  Utils::Vector3d
  get_node_velocity_at_boundary(const Utils::Vector3i &node) const {
    return m_callback.get_node_boundary_velocity(node);
  }

  void set_node_velocity_at_boundary(Utils::Vector3i const &node,
                                     Utils::Vector3d const &v,
                                     BlockAndCell const &bc) {
    auto [flag_field, boundary_flag] = get_flag_field_and_flag(bc.block);
    m_callback.set_node_boundary_velocity(node, v);
    flag_field->addFlag(bc.cell, boundary_flag);
  }

  void remove_node_from_boundary(Utils::Vector3i const &node,
                                 BlockAndCell const &bc) {
    auto [flag_field, boundary_flag] = get_flag_field_and_flag(bc.block);
    m_callback.unset_node_boundary_velocity(node);
    flag_field->removeFlag(bc.cell, boundary_flag);
  }

  /** Assign velocity boundary conditions to boundary cells. */
  void ubb_update() {
    m_boundary->template fillFromFlagField<FlagField>(
        m_blocks, m_flag_field_id, Boundary_flag, Domain_flag);
  }

private:
  std::shared_ptr<StructuredBlockForest> m_blocks;
  BlockDataID m_pdf_field_id;
  BlockDataID m_flag_field_id;
  DynamicVelocityCallback m_callback;
  std::shared_ptr<Dynamic_UBB> m_boundary;

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
