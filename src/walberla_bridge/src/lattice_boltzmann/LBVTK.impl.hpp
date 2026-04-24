/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

/**
 * @file
 * Out-of-class VTK writer registration definition for
 * @ref walberla::LBWalberlaImpl.
 */

#include <field/iterators/IteratorMacros.h>

#include <memory>
#include <optional>
#include <string>

namespace walberla {

/**
 * @brief Base class for LB field VTK writers.
 * Provides unit conversion and field access for cell-based VTK output.
 * All field is copied to a Vector container before writing.
 * @tparam FloatType   Internal LB precision (float or double).
 * @tparam VecType     Vector type to copy the field data into.
 * @tparam F_SIZE_ARG  Number of components per cell (1, 3, or 9).
 * @tparam OutputType  VTK output precision (default: float).
 */
template <typename FloatType, typename VecType, uint_t F_SIZE_ARG,
          typename OutputType>
class VTKWriter : public vtk::BlockCellDataWriter<OutputType, F_SIZE_ARG> {
public:
  VTKWriter(ConstBlockDataID const &block_id, std::string const &id,
            FloatType unit_conversion)
      : vtk::BlockCellDataWriter<OutputType, F_SIZE_ARG>(id),
        m_conversion(unit_conversion), m_content{} {}

protected:
  void configure() override { WALBERLA_ASSERT_NOT_NULLPTR(this->block_); }

  std::size_t get_first_index(cell_idx_t const x, cell_idx_t const y,
                              cell_idx_t const z) {
    return (static_cast<std::size_t>(x) * m_dims[2] * m_dims[1] +
            static_cast<std::size_t>(y) * m_dims[2] +
            static_cast<std::size_t>(z)) *
           F_SIZE_ARG;
  }

  FloatType m_conversion;
  VecType m_content;
  Vector3<uint_t> m_dims;

public:
  void set_content(VecType content) { m_content = std::move(content); }

  void set_dims(Vector3<uint_t> dims) { m_dims = dims; }
};

template <typename FloatType, typename OutputType = float>
class DensityVTKWriter
    : public VTKWriter<FloatType, std::vector<FloatType>, 1u, OutputType> {
public:
  using Base = VTKWriter<FloatType, std::vector<FloatType>, 1u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const) override {
    WALBERLA_ASSERT(!this->m_content.empty());
    auto const density = this->m_content[this->get_first_index(x, y, z)];
    return numeric_cast<OutputType>(this->m_conversion * density);
  }
};

template <typename FloatType, typename OutputType = float>
class VelocityVTKWriter
    : public VTKWriter<FloatType, std::vector<FloatType>, 3u, OutputType> {
public:
  using Base = VTKWriter<FloatType, std::vector<FloatType>, 3u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const f) override {
    WALBERLA_ASSERT(!this->m_content.empty());
    auto velocity = this->m_content[this->get_first_index(x, y, z) + f];
    return numeric_cast<OutputType>(this->m_conversion * velocity);
  }
};

template <typename FloatType, typename OutputType = float>
class PressureTensorVTKWriter
    : public VTKWriter<FloatType, std::vector<FloatType>, 9u, OutputType> {
public:
  using Base = VTKWriter<FloatType, std::vector<FloatType>, 9u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const f) override {
    WALBERLA_ASSERT(!this->m_content.empty());
    auto pressure = this->m_content[this->get_first_index(x, y, z) + f];
    return numeric_cast<OutputType>(this->m_conversion * pressure);
  }
};

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::register_vtk_field_writers(
    walberla::vtk::VTKOutput &vtk_obj, LatticeModel::units_map const &units,
    int flag_observables) {
  if (flag_observables & static_cast<int>(OutputVTK::density)) {
    auto const unit_conversion = FloatType_c(units.at("density"));
    auto const blocks = m_lattice->get_blocks();
    WALBERLA_ASSERT_NOT_NULLPTR(blocks);
    auto density_writer = std::make_shared<DensityVTKWriter<FloatType, float>>(
        m_pdf_field_id, "density", unit_conversion);
    auto before_function = [this, blocks, density_writer]() {
      for (auto &block : *blocks) {
        auto *pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        auto const bci = pdf_field->xyzSize();
        density_writer->set_content(
            lbm::accessor::Density::get(pdf_field, m_density, bci));
        density_writer->set_dims(
            Vector3<uint_t>(bci.xSize(), bci.ySize(), bci.zSize()));
      }
    };
    vtk_obj.addBeforeFunction(std::move(before_function));
    vtk_obj.addCellDataWriter(density_writer);
  }
  if (flag_observables & static_cast<int>(OutputVTK::velocity_vector)) {
    auto const unit_conversion = FloatType_c(units.at("velocity"));
    auto const blocks = m_lattice->get_blocks();
    WALBERLA_ASSERT_NOT_NULLPTR(blocks);
    auto velocity_writer =
        std::make_shared<VelocityVTKWriter<FloatType, float>>(
            m_pdf_field_id, "velocity_vector", unit_conversion);
    auto before_function = [this, blocks, velocity_writer]() {
      for (auto &block : *blocks) {
        auto *velocity_field =
            block.template getData<VectorField>(m_velocity_field_id);
        auto const bci = velocity_field->xyzSize();
        velocity_writer->set_content(
            lbm::accessor::Vector::get(velocity_field, bci));
        velocity_writer->set_dims(
            Vector3<uint_t>(bci.xSize(), bci.ySize(), bci.zSize()));
      }
    };
    vtk_obj.addBeforeFunction(std::move(before_function));
    vtk_obj.addCellDataWriter(velocity_writer);
  }
  if (flag_observables & static_cast<int>(OutputVTK::pressure_tensor)) {
    auto const unit_conversion = FloatType_c(units.at("pressure"));
    auto const blocks = m_lattice->get_blocks();
    WALBERLA_ASSERT_NOT_NULLPTR(blocks);
    auto pressure_writer =
        std::make_shared<PressureTensorVTKWriter<FloatType, float>>(
            m_pdf_field_id, "pressure_tensor", unit_conversion);
    auto before_function = [this, blocks, pressure_writer]() {
      for (auto &block : *blocks) {
        auto *pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        auto const bci = pdf_field->xyzSize();
        auto values =
            lbm::accessor::PressureTensor::get(pdf_field, m_density, bci);
        for (std::size_t n = 0u; n < values.size(); n += 9u) {
          pressure_tensor_correction(
              std::span<FloatType, 9ul>(&values[n], 9ul));
        }
        pressure_writer->set_content(std::move(values));
        pressure_writer->set_dims(
            Vector3<uint_t>(bci.xSize(), bci.ySize(), bci.zSize()));
      }
    };
    vtk_obj.addBeforeFunction(std::move(before_function));
    vtk_obj.addCellDataWriter(pressure_writer);
  }
}

} // namespace walberla
