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

/**
 * @file
 * @ref walberla::LBWalberlaImpl implements the interface of the LB
 * waLBerla bridge. It is a templated class that is specialized by lattice
 * models created by lbmpy (see <tt>maintainer/walberla_kernels</tt>).
 */

#ifndef WALBERLA_BRIDGE_RESET_FORCE_HPP
#define WALBERLA_BRIDGE_RESET_FORCE_HPP

#include "blockforest/Initialization.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "boundary/BoundaryHandling.h"
#include "field/GhostLayerField.h"
#include "field/adaptors/GhostLayerFieldAdaptor.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/vtk/all.h"
#include "timeloop/SweepTimeloop.h"

#include "domain_decomposition/SharedSweep.h"

#include "lbm/sweeps/CellwiseSweep.h"

namespace walberla {

/** Sweep that swaps force_to_be_applied and last_applied_force
and resets force_to_be_applied to the global external force
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
    PdfField *pdf_field = block->template getData<PdfField>(m_pdf_field_id);
    ForceField *force_field =
        block->template getData<ForceField>(m_last_applied_force_field_id);
    ForceField *force_to_be_applied =
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
