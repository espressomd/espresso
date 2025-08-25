/*
 * Copyright (C) 2025 The ESPResSo project
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
#include "../../../../src/electrokinetics/generated_kernels/EK_FieldAccessors_double_precision_CUDA.cuh"
#include "../../../../src/electrokinetics/generated_kernels/EK_FieldAccessors_single_precision.h"
#include "../../../../src/electrokinetics/generated_kernels/EK_FieldAccessors_single_precision_CUDA.cuh"
#include "../../BlockAndCell.hpp"

#include <domain_decomposition/BlockDataID.h>
#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>
#include <field/vtk/VTKWriter.h>
#include <gpu/HostFieldAllocator.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <utility>

namespace walberla {
template <typename FloatType> class FFT_CUDA;

template <typename FloatType> class FFT_GPU : public PoissonSolver {
private:
  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  std::shared_ptr<FFT_CUDA<FloatType>> fft_cuda;
  using PotentialField = gpu::GPUField<FloatType>;
  using PotentialFieldCpu = GhostLayerField<FloatType, 1>;

  std::optional<BlockDataID> m_potential_cpu_field_id;
  std::shared_ptr<gpu::HostFieldAllocator<FloatType>> m_host_field_allocator;

public:
  FFT_GPU() = default;
  FFT_GPU(std::shared_ptr<LatticeWalberla> lattice, double permittivity)
      : PoissonSolver(lattice, permittivity) {
    fft_cuda = std::make_shared<FFT_CUDA<FloatType>>(lattice, permittivity);
    m_host_field_allocator =
        std::make_shared<gpu::HostFieldAllocator<FloatType>>();
  }
  ~FFT_GPU() override = default;

  void reset_charge_field() override { fft_cuda->reset_charge_field(); }

  void add_charge_to_field(std::size_t id, double valency,
                           bool is_double_precision) override {
    fft_cuda->add_charge_to_field(id, valency, is_double_precision);
  }

  std::size_t get_potential_field_id() const noexcept override {
    return fft_cuda->get_potential_field_id();
  }

  void solve() override {
    fft_cuda->solve();
    integrate_vtk_writers();
  }

  void set_permittivity(double permittivity) noexcept override {
    fft_cuda->set_permittivity(permittivity);
  }

  [[nodiscard]] double get_permittivity() const noexcept override {
    return fft_cuda->get_permittivity();
  }

  [[nodiscard]] LatticeWalberla const &get_lattice() const noexcept override {
    return fft_cuda->get_lattice();
  }

  [[nodiscard]] std::optional<double>
  get_node_potential(Utils::Vector3i const &node,
                     bool consider_ghosts = false) override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);

    if (!bc || (get_potential_field_id() == 0))
      return std::nullopt;

    auto const potential_field = bc->block->template getData<PotentialField>(
        domain_decomposition::BlockDataID(get_potential_field_id()));
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
          auto const potential_field = block.template getData<PotentialField>(
              domain_decomposition::BlockDataID(get_potential_field_id()));
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
  class PotentialVTKWriter
      : public VTKWriter<PotentialFieldCpu, 1u, OutputType> {
  public:
    using Base = VTKWriter<PotentialFieldCpu, 1u, OutputType>;
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
    auto const allocate_cpu_field_if_empty =
        [&]<typename Field>(auto const &blocks, std::string name,
                            std::optional<BlockDataID> &cpu_field) {
          if (not cpu_field) {
            cpu_field = field::addToStorage<Field>(
                blocks, name, FloatType{0}, field::fzyx,
                get_lattice().get_ghost_layers(), m_host_field_allocator);
          }
        };
    if (flag_observables & static_cast<int>(EKPoissonOutputVTK::potential)) {
      auto const unit_conversion = FloatType_c(units.at("potential"));
      auto const &blocks = get_lattice().get_blocks();
      allocate_cpu_field_if_empty.template operator()<PotentialFieldCpu>(
          blocks, "potential_cpu", m_potential_cpu_field_id);
      vtk_obj.addBeforeFunction(
          gpu::fieldCpyFunctor<PotentialFieldCpu, PotentialField>(
              blocks, *m_potential_cpu_field_id,
              domain_decomposition::BlockDataID(get_potential_field_id())));
      vtk_obj.addCellDataWriter(make_shared<PotentialVTKWriter<float>>(
          *m_potential_cpu_field_id, "potential", unit_conversion));
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
  void ghost_communication() { fft_cuda->ghost_communication(); }
};

} // namespace walberla
