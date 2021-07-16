/*
 * Copyright (C) 2020 The ESPResSo project
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

#ifndef WALBERLA_BRIDGE_RESET_FORCE_HPP
#define WALBERLA_BRIDGE_RESET_FORCE_HPP

#include "domain_decomposition/SharedSweep.h"

#include "lbm/sweeps/CellwiseSweep.h"

#include "walberla_utils.hpp"

#include <utils/Vector.hpp>

namespace walberla {

/** Sweep that swaps @c force_to_be_applied and @c last_applied_force
 *  and resets @c force_to_be_applied to the global external force.
 */
template <typename PdfField, typename ForceField> class ResetForce {
public:
  ResetForce(const BlockDataID &pdf_field_id,
             const BlockDataID &last_applied_force_field_id,
             const BlockDataID &force_to_be_applied_id)
      : m_pdf_field_id(pdf_field_id),
        m_last_applied_force_field_id(last_applied_force_field_id),
        m_force_to_be_applied_id(force_to_be_applied_id),
        m_ext_force(Vector3<real_t>{0, 0, 0}){};

  void set_ext_force(const Utils::Vector3d &ext_force) {
    m_ext_force = to_vector3(ext_force);
  }

  Utils::Vector3d get_ext_force() const { return to_vector3d(m_ext_force); };

  void operator()(IBlock *block) {
    auto *pdf_field = block->template getData<PdfField>(m_pdf_field_id);
    auto *force_field =
        block->template getData<ForceField>(m_last_applied_force_field_id);
    auto *force_to_be_applied =
        block->template getData<ForceField>(m_force_to_be_applied_id);

    force_field->swapDataPointers(force_to_be_applied);

    WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(force_field, {
      Cell cell(x, y, z);
      for (int i : {0, 1, 2})
        force_field->get(x, y, z, i) += m_ext_force[i];
    });
    WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(force_to_be_applied, {
      Cell cell(x, y, z);
      for (int i : {0, 1, 2})
        force_to_be_applied->get(cell, i) = real_t{0};
    });
  }

private:
  const BlockDataID m_pdf_field_id;
  const BlockDataID m_last_applied_force_field_id;
  const BlockDataID m_force_to_be_applied_id;
  Vector3<real_t> m_ext_force;
};
} // namespace walberla
#endif
