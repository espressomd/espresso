/*
 * Copyright (C) 2020-2023 The ESPResSo project
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

#include "generated_kernels/FieldAccessorsDoublePrecision.h"
#include "generated_kernels/FieldAccessorsSinglePrecision.h"

#include <walberla_bridge/utils/walberla_utils.hpp>

#include <core/math/Vector3.h>
#include <domain_decomposition/SharedSweep.h>
#include <lbm/sweeps/CellwiseSweep.h>

#include <utils/Vector.hpp>

namespace walberla {

/** Sweep that swaps @c force_to_be_applied and @c last_applied_force
 *  and resets @c force_to_be_applied to the global external force.
 */
template <typename PdfField, typename ForceField> class ResetForce {
  using FloatType = typename PdfField::value_type;

public:
  ResetForce(BlockDataID const &last_applied_force_field_id,
             BlockDataID const &force_to_be_applied_id)
      : m_last_applied_force_field_id(last_applied_force_field_id),
        m_force_to_be_applied_id(force_to_be_applied_id),
        m_ext_force(Vector3<FloatType>{0, 0, 0}) {}

  void set_ext_force(Utils::Vector3d const &ext_force) {
    m_ext_force = to_vector3<FloatType>(ext_force);
  }

  Utils::Vector3d get_ext_force() const { return to_vector3d(m_ext_force); }

  void operator()(IBlock *block) {
    auto force_field =
        block->template getData<ForceField>(m_last_applied_force_field_id);
    auto force_to_be_applied =
        block->template getData<ForceField>(m_force_to_be_applied_id);

    force_field->swapDataPointers(force_to_be_applied);

    lbm::accessor::Vector::add_to_all(force_field, m_ext_force);
    lbm::accessor::Vector::broadcast(force_to_be_applied,
                                     Vector3<FloatType>{0});
  }

private:
  const BlockDataID m_last_applied_force_field_id;
  const BlockDataID m_force_to_be_applied_id;
  Vector3<FloatType> m_ext_force;
};

} // namespace walberla
