/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "PoissonSolver.hpp"

#include "../../../../src/electrokinetics/generated_kernels/EK_FieldAccessors_double_precision.h"
#include "../../../../src/electrokinetics/generated_kernels/EK_FieldAccessors_single_precision.h"
#include "../../BlockAndCell.hpp"

#include <blockforest/communication/UniformBufferedScheme.h>
#include <domain_decomposition/BlockDataID.h>
#include <fft/Fft.h>
#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>
#include <field/communication/PackInfo.h>
#include <field/vtk/VTKWriter.h>
#include <stencil/D3Q27.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <utility>

namespace walberla {

template <typename FloatType> class PoissonSolverFFT : public PoissonSolver {
private:
  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  domain_decomposition::BlockDataID m_potential_field_id;

  using PotentialField = GhostLayerField<FloatType, 1>;

  std::shared_ptr<fft::FourierTransform<PotentialField>> m_ft;
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

  using FullCommunicator =
      blockforest::communication::UniformBufferedScheme<stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;

public:
  PoissonSolverFFT(std::shared_ptr<LatticeWalberla> lattice,
                   double permittivity)
      : PoissonSolver(std::move(lattice), permittivity) {
    m_blocks = get_lattice().get_blocks();

    Vector3<uint_t> dim(m_blocks->getNumberOfXCells(),
                        m_blocks->getNumberOfYCells(),
                        m_blocks->getNumberOfZCells());
    auto const greens = [dim](uint_t x, uint_t y, uint_t z) -> real_t {
      if (x == 0u && y == 0u && z == 0u)
        return 0.;
      return -0.5 /
             (std::cos(2. * std::numbers::pi * real_c(x) / real_c(dim[0])) +
              std::cos(2. * std::numbers::pi * real_c(y) / real_c(dim[1])) +
              std::cos(2. * std::numbers::pi * real_c(z) / real_c(dim[2])) -
              3.) /
             real_c(dim[0] * dim[1] * dim[2]);
    };

    m_potential_field_id = field::addToStorage<PotentialField>(
        get_lattice().get_blocks(), "potential field", 0.0, field::fzyx,
        get_lattice().get_ghost_layers());

    m_ft = std::make_shared<fft::FourierTransform<PotentialField>>(
        m_blocks, m_potential_field_id, greens);

    m_full_communication =
        std::make_shared<FullCommunicator>(get_lattice().get_blocks());
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PotentialField>>(
            m_potential_field_id));
  }
  ~PoissonSolverFFT() override = default;

  void reset_charge_field() override {
    // the FFT-solver re-uses the potential field for the charge
    auto const potential_id = walberla::BlockDataID(get_potential_field_id());

    for (auto &block : *get_lattice().get_blocks()) {
      auto field = block.template getData<PotentialField>(potential_id);
      WALBERLA_FOR_ALL_CELLS_XYZ(field, field->get(x, y, z) = 0.;)
    }
  }

  void add_charge_to_field(std::size_t id, double valency,
                           bool is_double_precision) override {
    auto const factor = FloatType_c(valency) / FloatType_c(get_permittivity());
    // the FFT-solver re-uses the potential field for the charge
    const auto charge_id = walberla::BlockDataID(get_potential_field_id());
    const auto density_id = walberla::BlockDataID(id);
    for (auto &block : *get_lattice().get_blocks()) {
      auto charge_field = block.template getData<PotentialField>(charge_id);
      if (is_double_precision) {
        auto density_field =
            block.template getData<walberla::GhostLayerField<double, 1>>(
                density_id);
        WALBERLA_FOR_ALL_CELLS_XYZ(
            charge_field, charge_field->get(x, y, z) +=
                          factor * FloatType_c(density_field->get(x, y, z));)
      } else {
        auto density_field =
            block.template getData<walberla::GhostLayerField<float, 1>>(
                density_id);
        WALBERLA_FOR_ALL_CELLS_XYZ(
            charge_field, charge_field->get(x, y, z) +=
                          factor * FloatType_c(density_field->get(x, y, z));)
      }
    }
  }

  [[nodiscard]] std::size_t get_potential_field_id() const noexcept override {
    return static_cast<std::size_t>(m_potential_field_id);
  }

  void solve() override {
    (*m_ft)();
    ghost_communication();
    integrate_vtk_writers();
  }

  [[nodiscard]] std::optional<double>
  get_node_potential(Utils::Vector3i const &node,
                     bool consider_ghosts = false) override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);

    if (!bc || (get_potential_field_id() == 0))
      return std::nullopt;

    auto const potential_field =
        bc->block->template getData<PotentialField>(m_potential_field_id);
    return {double_c(
        walberla::ek::accessor::Scalar::get(potential_field, bc->cell))};
  }

  [[nodiscard]] std::vector<double>
  get_slice_potential(Utils::Vector3i const &lower_corner,
                      Utils::Vector3i const &upper_corner) const override {
    std::vector<double> out;
    uint_t values_size = 0;
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      out = std::vector<double>(ci->numCells());
      for (auto &block : *lattice.get_blocks()) {
        auto const block_offset = lattice.get_block_corner(block, true);
        if (auto const bci = get_block_interval(
                lattice, lower_corner, upper_corner, block_offset, block)) {
          auto const potential_field =
              block.template getData<PotentialField>(m_potential_field_id);
          auto const values = ek::accessor::Scalar::get(potential_field, *bci);
          assert(values.size() == bci->numCells());
          values_size += bci->numCells();
          auto kernel = [&values, &out](unsigned const block_index,
                                        unsigned const local_index,
                                        Utils::Vector3i const &node) {
            out[local_index] = double_c(values[block_index]);
          };

          copy_block_buffer(*bci, *ci, block_offset, lower_corner, kernel);
        }
      }
      assert(values_size == ci->numCells());
    }
    return out;
  }

protected:
  template <typename Field_T, uint_t F_SIZE_ARG, typename OutputType>
  class VTKWriter : public vtk::BlockCellDataWriter<OutputType, F_SIZE_ARG> {
  public:
    VTKWriter(ConstBlockDataID const &block_id, std::string const &id,
              FloatType unit_conversion)
        : vtk::BlockCellDataWriter<OutputType, F_SIZE_ARG>(id),
          m_block_id(block_id), m_field(nullptr),
          m_conversion(unit_conversion) {}

  protected:
    void configure() override {
      WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
      m_field = this->block_->template getData<Field_T>(m_block_id);
    }

    ConstBlockDataID const m_block_id;
    Field_T const *m_field;
    FloatType const m_conversion;
  };

  template <typename OutputType = float>
  class PotentialVTKWriter : public VTKWriter<PotentialField, 1u, OutputType> {
  public:
    using Base = VTKWriter<PotentialField, 1u, OutputType>;
    using Base::Base;
    using Base::evaluate;

  protected:
    OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                        cell_idx_t const z, cell_idx_t const) override {
      WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
      auto const potential =
          walberla::ek::accessor::Scalar::get(this->m_field, {x, y, z});
      return numeric_cast<OutputType>(this->m_conversion * potential);
    }
  };

  void register_vtk_field_writers(walberla::vtk::VTKOutput &vtk_obj,
                                  LatticeModel::units_map const &units,
                                  int flag_observables) override {
    if (flag_observables & static_cast<int>(EKPoissonOutputVTK::potential)) {
      auto const unit_conversion = FloatType_c(units.at("potential"));
      vtk_obj.addCellDataWriter(make_shared<PotentialVTKWriter<float>>(
          m_potential_field_id, "potential", unit_conversion));
    }
  }

  void integrate_vtk_writers() override {
    for (auto const &it : m_vtk_auto) {
      auto &vtk_handle = it.second;
      if (vtk_handle->enabled) {
        vtk::writeFiles(vtk_handle->ptr)();
        vtk_handle->execution_count++;
      }
    }
  }

private:
  void ghost_communication() { (*m_full_communication)(); }
};

} // namespace walberla
