/*
 * Copyright (C) 2020-2021 The ESPResSo project
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

#include "generated_kernels/macroscopic_values_accessors_double_precision.h"
#include "generated_kernels/macroscopic_values_accessors_single_precision.h"

namespace walberla {
namespace lbm {

// code adapted from "lbm/vtk/Density.h"
template <typename LatticeModel_T, typename OutputType = float>
class DensityVTKWriter : public vtk::BlockCellDataWriter<OutputType, 1> {
public:
  using PdfField_T = typename LatticeModel_T::PdfField;

  DensityVTKWriter(ConstBlockDataID const &pdf, std::string const &id)
      : vtk::BlockCellDataWriter<OutputType, 1>(id), bdid_(pdf), pdf_(nullptr) {
  }

protected:
  void configure() override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
    pdf_ = this->block_->template getData<PdfField_T>(bdid_);
  }

  OutputType evaluate(const cell_idx_t x, const cell_idx_t y,
                      const cell_idx_t z, const cell_idx_t /*f*/) override {
    WALBERLA_ASSERT_NOT_NULLPTR(pdf_);
    return numeric_cast<OutputType>(
        lbm::accessor::Density::get(*pdf_, x, y, z));
  }

  ConstBlockDataID const bdid_;
  PdfField_T const *pdf_;
};

// code adapted from "lbm/vtk/PressureTensor.h"
template <typename LatticeModel_T, typename OutputType = float>
class PressureTensorVTKWriter : public vtk::BlockCellDataWriter<OutputType, 9> {
public:
  using PdfField_T = typename LatticeModel_T::PdfField;
  using FloatType = typename PdfField_T::value_type;

  PressureTensorVTKWriter(ConstBlockDataID const &pdf, std::string const &id,
                          FloatType off_diag_factor)
      : vtk::BlockCellDataWriter<OutputType, 9>(id), bdid_(pdf), pdf_(nullptr),
        off_diag_factor_(off_diag_factor) {}

protected:
  void configure() override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
    pdf_ = this->block_->template getData<PdfField_T>(bdid_);
  }

  OutputType evaluate(const cell_idx_t x, const cell_idx_t y,
                      const cell_idx_t z, const cell_idx_t f) override {
    WALBERLA_ASSERT_NOT_NULLPTR(pdf_);
    Matrix3<typename PdfField_T::value_type> pressureTensor;
    lbm::accessor::PressureTensor::get(pressureTensor, *pdf_, x, y, z);
    auto const revert_factor =
        (f == 0 or f == 4 or f == 8) ? FloatType{1} : off_diag_factor_;
    return numeric_cast<OutputType>(revert_factor * pressureTensor[f]);
  }

  ConstBlockDataID const bdid_;
  PdfField_T const *pdf_;
  FloatType const off_diag_factor_;
};

} // namespace lbm
} // namespace walberla
