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

#include "vtk/BlockCellDataWriter.h"

#include "generated_kernels/macroscopic_values_accessors.h"

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

  PressureTensorVTKWriter(ConstBlockDataID const &pdf, std::string const &id)
      : vtk::BlockCellDataWriter<OutputType, 9>(id), bdid_(pdf), pdf_(nullptr) {
  }

protected:
  void configure() override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
    pdf_ = this->block_->template getData<PdfField_T>(bdid_);
  }

  OutputType evaluate(const cell_idx_t x, const cell_idx_t y,
                      const cell_idx_t z, const cell_idx_t f) override {
    WALBERLA_ASSERT_NOT_NULLPTR(pdf_);
    Matrix3<real_t> pressureTensor;
    lbm::accessor::PressureTensor::get(pressureTensor, *pdf_, x, y, z);
    return numeric_cast<OutputType>(pressureTensor[f]);
  }

  ConstBlockDataID const bdid_;
  PdfField_T const *pdf_;
};

} // namespace lbm
} // namespace walberla
