/*
 * Copyright (C) 2020-2023 The ESPResSo project
 * Copyright (C) 2017-2021 The waLBerla project
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

#include <vtk/BlockCellDataWriter.h>

#include "../generated_kernels/FieldAccessorsDoublePrecision.h"
#include "../generated_kernels/FieldAccessorsSinglePrecision.h"

namespace walberla {
namespace lbm {

template <typename LatticeModel_T, typename OutputType = float>
class VelocityVTKWriter : public vtk::BlockCellDataWriter<OutputType, 3> {
public:
  using Field_T = typename LatticeModel_T::VectorField;
  using FloatType = typename Field_T::value_type;

  VelocityVTKWriter(ConstBlockDataID const &vel_id, std::string const &id,
                    FloatType unit_conversion)
      : vtk::BlockCellDataWriter<OutputType, 3>(id), m_block_id(vel_id),
        m_vel_field(nullptr), m_conversion(unit_conversion) {}

protected:
  void configure() override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
    m_vel_field = this->block_->template getData<Field_T>(m_block_id);
  }

  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const f) override {
    WALBERLA_ASSERT_NOT_NULLPTR(m_vel_field);
    auto const velocity = lbm::accessor::Vector::get(m_vel_field, {x, y, z});
    return numeric_cast<OutputType>(m_conversion * velocity[f]);
  }

  ConstBlockDataID const m_block_id;
  Field_T const *m_vel_field;
  FloatType const m_conversion;
};

} // namespace lbm
} // namespace walberla
