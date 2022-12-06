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

#include "LeesEdwardsPack.hpp"

#include "blockforest/StructuredBlockForest.h"
#include "domain_decomposition/all.h"
#include "stencil/D3Q19.h"

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>

namespace walberla {

/** Lees-Edwards sweep.
 * TODO WALBERLA
 *   Currently only works for 1 MPI rank! It should work in parallel if the MPI
 *   domain decomposition for the structured block forest doesn't partition
 *   along the shear direction. For example if the shear direction goes along
 *   the z-axis, it should be possible to run on 4 MPI ranks with [2, 2, 1].
 *   At the moment, ESPResSo requires system.cell_system.node_grid to be in
 *   decreasing order, therefore parallelization requires a shear direction
 *   along the z-axis and a MPI node_grid of [x, y, 1] with x >= y. This
 *   restriction on the ordering of the node_grid may be lifted in the
 *   distant future, when our FFT algorithm is replaced by a new one.
 */
template <class FieldType, typename FloatType>
class InterpolateAndShiftAtBoundary {
public:
  InterpolateAndShiftAtBoundary(
      std::shared_ptr<StructuredBlockForest> blocks, BlockDataID field_id,
      BlockDataID tmp_field_id, unsigned int n_ghost_layers,
      unsigned int shear_direction, unsigned int shear_plane_normal,
      std::function<double()> get_pos_offset,
      std::function<double()> get_shift = []() { return 0.0; })
      : m_blocks(std::move(blocks)), m_field_id(field_id),
        m_tmp_field_id(tmp_field_id), m_n_ghost_layers(uint_c(n_ghost_layers)),
        m_shear_direction(uint_c(shear_direction)),
        m_shear_plane_normal(uint_c(shear_plane_normal)),
        m_get_pos_offset(std::move(get_pos_offset)),
        m_get_shift(std::move(get_shift)) {
    if (m_n_ghost_layers != 1u) {
      throw std::runtime_error("The Lees-Edwards sweep is implemented "
                               "for a ghost layer of thickness 1");
    }
    if (m_shear_plane_normal == 0u) {
      m_slab_min = stencil::W;
      m_slab_max = stencil::E;
    } else if (m_shear_plane_normal == 1u) {
      m_slab_min = stencil::S;
      m_slab_max = stencil::N;
    } else if (m_shear_plane_normal == 2u) {
      m_slab_min = stencil::B;
      m_slab_max = stencil::T;
    }
  }
  typedef stencil::D3Q19 Stencil;

  FloatType get_pos_offset() const {
    return numeric_cast<FloatType>(m_get_pos_offset());
  }

  FloatType get_shift() const { return numeric_cast<FloatType>(m_get_shift()); }

  void operator()(IBlock *block) {
    kernel(block, m_slab_min);
    kernel(block, m_slab_max);
  }

private:
  void kernel(IBlock *block, stencil::Direction slab_dir) {
    // setup lengths
    assert(m_blocks->getNumberOfCells(*block, m_shear_plane_normal) >= 2u);
    auto const dir = m_shear_direction;
    auto const dim = cell_idx_c(m_blocks->getNumberOfCells(*block, dir));
    auto const length = numeric_cast<FloatType>(dim);
    auto const weight =
        std::abs(std::fmod(get_pos_offset() + length, FloatType{1}));

    // setup slab
    auto field = block->template getData<FieldType>(m_field_id);
    auto tmp_field = block->template getData<FieldType>(m_tmp_field_id);

    CellInterval ci;
    field->getGhostRegion(slab_dir, ci, cell_idx_t{1}, true);

    // shift
    auto const shift = get_shift();
    // Note that the offset is applied to the interpolation source rather than
    // the target
    auto const prefactor =
        ((slab_dir == m_slab_max) ? FloatType{-1} : FloatType{1});
    auto const offset = get_pos_offset() * prefactor;
    for (auto const &&cell : ci) {
      Cell source1 = cell;
      Cell source2 = cell;
      source1[dir] = cell_idx_c(std::floor(
                         static_cast<FloatType>(source1[dir]) + offset)) %
                     dim;
      source1[dir] = cell_idx_c(static_cast<FloatType>(source1[dir]) + length);
      source1[dir] = cell_idx_c(source1[dir] % dim);

      source2[dir] =
          cell_idx_c(std::ceil(static_cast<FloatType>(source2[dir]) + offset)) %
          dim;
      source2[dir] = cell_idx_c(static_cast<FloatType>(source2[dir]) + length);
      source2[dir] = cell_idx_c(source2[dir] % dim);

      for (uint_t q = 0; q < FieldType::F_SIZE; ++q) {
        tmp_field->get(cell, q) = field->get(source1, q) * (1 - weight) +
                                  field->get(source2, q) * weight;
      }
      tmp_field->get(cell, m_shear_direction) -= prefactor * shift;
    }

    // swap
    for (auto const &&cell : ci) {
      for (uint_t q = 0; q < FieldType::F_SIZE; ++q) {
        field->get(cell, q) = tmp_field->get(cell, q);
      }
    }
  }

private:
  std::shared_ptr<StructuredBlockForest> m_blocks;
  BlockDataID m_field_id;
  BlockDataID m_tmp_field_id;
  uint_t m_n_ghost_layers;
  uint_t m_shear_direction;
  uint_t m_shear_plane_normal;
  std::function<double()> m_get_pos_offset;
  std::function<double()> m_get_shift;
  stencil::Direction m_slab_min;
  stencil::Direction m_slab_max;
};
} // namespace walberla
