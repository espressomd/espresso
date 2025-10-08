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

#include "utils/types_conversion.hpp"

#include <blockforest/StructuredBlockForest.h>
#include <domain_decomposition/BlockDataID.h>
#include <domain_decomposition/IBlock.h>
#include <field/FlagField.h>
#include <field/FlagUID.h>

#include <utils/Vector.hpp>

#if defined(__CUDACC__)
#include <thrust/device_vector.h>
#endif

#include <cassert>
#include <functional>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace walberla {

/**
 * @brief Boundary class optimized for sparse data.
 *
 * Instead of storing the boundary data on a vector field,
 * store individual vectors in a map.
 * The global cell is used as key.
 *
 * Requires a custom communicator:
 * @ref walberla::field::communication::BoundaryPackInfo.
 */
template <typename FloatType, typename ValueType, typename BoundaryClass>
class BoundaryHandling {
private:
  /** Flag for domain cells, i.e. all cells. */
  FlagUID const Domain_flag{"domain"};
  /** Flag for boundary cells. */
  FlagUID const Boundary_flag{"boundary"};

  /** Container for the map between cells and values. */
  class DynamicValueCallback {
  public:
    DynamicValueCallback() {
      m_value_boundary =
          std::make_shared<typename decltype(m_value_boundary)::element_type>();
    }

    [[nodiscard]] ValueType operator()(
        Cell const &local,
        std::shared_ptr<blockforest::StructuredBlockForest> const &blocks,
        IBlock &block) const {
      Cell global;
      blocks->transformBlockLocalToGlobalCell(global, block, local);
      return get_value(global);
    }

    void set_node_boundary_value(Cell const &global, ValueType const &val) {
      (*m_value_boundary)[global] = val;
    }

    void unset_node_boundary_value(Cell const &global) {
      assert(m_value_boundary->contains(global));
      m_value_boundary->erase(global);
    }

    [[nodiscard]] auto &get_node_boundary_value(Cell const &global) const {
      return get_value(global);
    }

    bool node_is_boundary(Cell const &global) const {
      return m_value_boundary->contains(global);
    }

#if defined(__CUDACC__)
    /**
     * @brief Build a flattened version of the unordered map container.
     * The coordinate list (COO) format is used, similar to how sparse
     * matrices are compressed, although here zero is a valid value.
     * Indices are relative to the origin of the local halo.
     */
    void rebuild_flat_map_device(CellInterval const &local_domain) {
      std::vector<int> indices;
      std::vector<FloatType> values;
      auto const &local_origin = local_domain.min();
      for (auto const &[cell, value] : *m_value_boundary) {
        if (local_domain.contains(cell)) {
          for (auto i : {0, 1, 2}) {
            indices.emplace_back(cell[i] - local_origin[i]);
          }
          if constexpr (std::is_arithmetic_v<ValueType>) {
            values.emplace_back(static_cast<FloatType>(value));
          } else {
            for (auto i : {0, 1, 2}) {
              values.emplace_back(static_cast<FloatType>(value[i]));
            }
          }
        }
      }
      m_flat_indices = decltype(m_flat_indices)(indices.begin(), indices.end());
      m_flat_values = decltype(m_flat_values)(values.begin(), values.end());
    }

    auto get_flattened_map_device() const {
      return std::make_pair(&m_flat_indices, &m_flat_values);
    }
#endif

  private:
#if defined(__CUDACC__)
    thrust::device_vector<int> m_flat_indices;
    thrust::device_vector<FloatType> m_flat_values;
#endif
    std::shared_ptr<std::unordered_map<Cell, ValueType>> m_value_boundary;
    static constexpr ValueType default_value{};

    [[nodiscard]] auto const &get_value(Cell const &cell) const {
      if (m_value_boundary->contains(cell)) {
        return m_value_boundary->at(cell);
      }
      return default_value;
    }
  };

  [[nodiscard]] inline auto get_flag_field_and_flag(IBlock *block) const {
    auto const flag_field =
        block->template uncheckedFastGetData<FlagField>(m_flag_field_id);
    auto const boundary_flag = flag_field->getFlag(Boundary_flag);
    return std::make_tuple(flag_field, boundary_flag);
  }

public:
  using value_type = ValueType;
  using FlagField = field::FlagField<uint8_t>;

  BoundaryHandling(std::shared_ptr<StructuredBlockForest> blocks,
                   BlockDataID value_field_id, BlockDataID flag_field_id,
                   CellInterval const &local_domain)
      : m_blocks(std::move(blocks)), m_flag_field_id(flag_field_id),
        m_callback(DynamicValueCallback()), m_local_domain(local_domain),
        m_pending_changes(false) {
    // reinitialize the flag field
    for (auto &block : *m_blocks) {
      flag_reset_kernel(block.template getData<FlagField>(m_flag_field_id));
    }
    // instantiate the boundary sweep
    std::function callback = m_callback;
    m_boundary =
        std::make_shared<BoundaryClass>(m_blocks, value_field_id, callback);
  }

  void operator()(IBlock *block) { (*m_boundary)(block); }

  [[nodiscard]] bool
  node_is_boundary(signed_integral_vector auto const &node) const {
    return m_callback.node_is_boundary(to_cell(node));
  }

  [[nodiscard]] auto &
  get_node_value_at_boundary(signed_integral_vector auto const &node) const {
    return m_callback.get_node_boundary_value(to_cell(node));
  }

  void set_node_value_at_boundary(signed_integral_vector auto const &node,
                                  ValueType const &v, BlockAndCell const &bc) {
    auto [flag_field, boundary_flag] = get_flag_field_and_flag(bc.block);
    m_callback.set_node_boundary_value(to_cell(node), v);
    flag_field->addFlag(bc.cell, boundary_flag);
    m_pending_changes = true;
  }

  void unpack_node(signed_integral_vector auto const &node,
                   ValueType const &v) {
    m_callback.set_node_boundary_value(to_cell(node), v);
  }

  void remove_node_from_boundary(signed_integral_vector auto const &node,
                                 BlockAndCell const &bc) {
    auto [flag_field, boundary_flag] = get_flag_field_and_flag(bc.block);
    m_callback.unset_node_boundary_value(to_cell(node));
    flag_field->removeFlag(bc.cell, boundary_flag);
    m_pending_changes = true;
  }

  /** Assign boundary conditions to boundary cells. */
  void boundary_update() {
    if (m_pending_changes) {
      m_boundary->template fillFromFlagField<FlagField>(
          m_blocks, m_flag_field_id, Boundary_flag, Domain_flag);
#if defined(__CUDACC__)
      m_callback.rebuild_flat_map_device(m_local_domain);
#endif
      m_pending_changes = false;
    }
  }

  std::tuple<int, int, int> block_dims(IBlock const &block) const {
    auto const field = block.template getData<FlagField>(m_flag_field_id);
    return {static_cast<int>(field->xSize()), static_cast<int>(field->ySize()),
            static_cast<int>(field->zSize())};
  }

#if defined(__CUDACC__)
  auto get_flattened_map_device() const {
    return m_callback.get_flattened_map_device();
  }
#endif

  Vector3<double> get_total_force(IBlock *block) const {
    return m_boundary->getForce(block);
  }

  auto const &get_force_vector(IBlock *block) const {
    using ForceVector = BoundaryClass::ForceVector;
    auto const force_vector_id = m_boundary->getForceVectorID();
    auto *forceVector = block->getData<ForceVector>(force_vector_id);
    forceVector->syncCPU();
    return m_boundary->getForceVector(block);
  }
  auto const &get_index_vector(IBlock const *block) const {
    return m_boundary->getIndexVector(block);
  }

private:
  std::shared_ptr<StructuredBlockForest> m_blocks;
  BlockDataID m_flag_field_id;
  DynamicValueCallback m_callback;
  std::shared_ptr<BoundaryClass> m_boundary;
  CellInterval m_local_domain;
  bool m_pending_changes;

  /** Register flags and reset all cells. */
  void flag_reset_kernel(FlagField *flag_field) {
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
