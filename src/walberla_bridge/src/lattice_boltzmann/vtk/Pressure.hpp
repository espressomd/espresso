/*
 * Copyright (C) 2021-2023 The ESPResSo project
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
class PressureTensorVTKWriter : public vtk::BlockCellDataWriter<OutputType, 9> {
public:
  using Field_T = typename LatticeModel_T::PdfField;
  using FloatType = typename Field_T::value_type;

  PressureTensorVTKWriter(ConstBlockDataID const &pdf, std::string const &id,
                          FloatType unit_conversion, FloatType off_diag_factor)
      : vtk::BlockCellDataWriter<OutputType, 9>(id), m_block_id(pdf),
        m_pdf_field(nullptr), m_conversion(unit_conversion),
        m_off_diag_factor(off_diag_factor) {}

protected:
  void configure() override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
    m_pdf_field = this->block_->template getData<Field_T>(m_block_id);
  }

  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const f) override {
    WALBERLA_ASSERT_NOT_NULLPTR(m_pdf_field);
    auto const pressure =
        lbm::accessor::PressureTensor::get(m_pdf_field, {x, y, z});
    auto const revert_factor =
        (f == 0 or f == 4 or f == 8) ? FloatType{1} : m_off_diag_factor;
    return numeric_cast<OutputType>(m_conversion * revert_factor * pressure[f]);
  }

  ConstBlockDataID const m_block_id;
  Field_T const *m_pdf_field;
  FloatType const m_conversion;
  FloatType const m_off_diag_factor;
};

} // namespace lbm
} // namespace walberla
