/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/BlockAndCell.hpp>
#include <walberla_bridge/electrokinetics/PoissonSolver.hpp>

#include "greens_function.hpp"

#include "generated_kernels/EK_FieldAccessors_double_precision.h"
#include "generated_kernels/EK_FieldAccessors_single_precision.h"
#if defined(__CUDACC__)
#include "generated_kernels/EK_FieldAccessors_double_precision_CUDA.cuh"
#include "generated_kernels/EK_FieldAccessors_single_precision_CUDA.cuh"
#endif

#include <utils/Vector.hpp>

#include <blockforest/communication/UniformBufferedScheme.h>
#include <domain_decomposition/BlockDataID.h>
#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>
#include <field/vtk/VTKWriter.h>
#include <stencil/D3Q27.h>
#include <waLBerlaDefinitions.h>
#if defined(__CUDACC__)
#include <gpu/AddGPUFieldToStorage.h>
#include <gpu/FieldAccessor.h>
#include <gpu/FieldIndexing.h>
#include <gpu/GPUField.h>
#include <gpu/HostFieldAllocator.h>
#include <gpu/Kernel.h>
#include <gpu/communication/UniformGPUScheme.h>
#endif

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-conversion"
#pragma clang diagnostic ignored "-Wimplicit-float-conversion"
#elif defined(__GNUC__) or defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__) or defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace walberla {
template <typename FloatType, lbmpy::Arch Architecture>
class PoissonSolverNone : public PoissonSolver {
private:
  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

protected:
  template <typename FT, lbmpy::Arch AT = lbmpy::Arch::CPU> struct FieldTrait {
    using PotentialField = field::GhostLayerField<FT, 1u>;
  };

#if defined(__CUDACC__)
  template <typename FT> struct FieldTrait<FT, lbmpy::Arch::GPU> {
    using PotentialField = gpu::GPUField<FT>;
  };
#endif

public:
  using PotentialField = FieldTrait<FloatType, Architecture>::PotentialField;

private:
  BlockDataID m_potential_field_id;

public:
  ~PoissonSolverNone() override = default;
  explicit PoissonSolverNone(std::shared_ptr<LatticeWalberla> lattice)
      : PoissonSolver(std::move(lattice), 0.0) {
    auto blocks = get_lattice().get_blocks();
#if defined(__CUDACC__)
    m_potential_field_id = gpu::addGPUFieldToStorage<PotentialField>(
        blocks, "potential field", 1u, field::fzyx,
        get_lattice().get_ghost_layers());
#else  // __CUDACC__
    m_potential_field_id = field::addToStorage<PotentialField>(
        blocks, "potential field", 0., field::fzyx,
        get_lattice().get_ghost_layers());
#endif // __CUDACC__
  }

  void setup_fft(bool use_gpu_aware) override {}

  [[nodiscard]] bool is_gpu() const noexcept override {
    return Architecture == lbmpy::Arch::GPU;
  }

  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same_v<FloatType, double>;
  }

  std::size_t get_potential_field_id() const noexcept override {
    return static_cast<std::size_t>(m_potential_field_id);
  }

  [[nodiscard]] std::optional<double>
  get_node_potential(Utils::Vector3i const &node,
                     bool consider_ghosts = false) override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);

    if (not bc or get_potential_field_id() == 0u)
      return std::nullopt;

    auto const potential_field =
        bc->block->template getData<PotentialField>(m_potential_field_id);
    return {double_c(
        walberla::ek::accessor::Scalar::get(potential_field, bc->cell))};
  }

  bool set_node_potential(Utils::Vector3i const &node,
                          double potential) override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc) {
      return false;
    }
    auto potential_field =
        bc->block->template getData<PotentialField>(m_potential_field_id);
    ek::accessor::Scalar::set(potential_field, FloatType_c(potential),
                              bc->cell);
    return true;
  }

  [[nodiscard]] std::vector<double>
  get_slice_potential(Utils::Vector3i const &lower_corner,
                      Utils::Vector3i const &upper_corner) const override {
    std::vector<double> out;
#ifndef NDEBUG
    uint_t values_size{0u};
#endif
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
#ifndef NDEBUG
          values_size += bci->numCells();
#endif
          auto kernel = [&values, &out](unsigned const block_index,
                                        unsigned const local_index,
                                        Utils::Vector3i const &) {
            out[local_index] = double_c(values[block_index]);
          };

          copy_block_buffer(*bci, *ci, block_offset, lower_corner, kernel);
        }
      }
      assert(values_size == ci->numCells());
    }
    return out;
  }

  void set_slice_potential(Utils::Vector3i const &lower_corner,
                           Utils::Vector3i const &upper_corner,
                           std::vector<double> const &potential) override {
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      assert(potential.size() == ci->numCells());
      for (auto &block : *lattice.get_blocks()) {
        auto const block_offset = lattice.get_block_corner(block, true);
        if (auto const bci = get_block_interval(
                lattice, lower_corner, upper_corner, block_offset, block)) {
          auto potential_field =
              block.template getData<PotentialField>(m_potential_field_id);
          std::vector<FloatType> values(bci->numCells());

          auto kernel = [&values, &potential](unsigned const block_index,
                                              unsigned const local_index,
                                              Utils::Vector3i const &) {
            values[block_index] =
                numeric_cast<FloatType>(potential[local_index]);
          };

          copy_block_buffer(*bci, *ci, block_offset, lower_corner, kernel);
          ek::accessor::Scalar::set(potential_field, values, *bci);
        }
      }
    }
  }
  void ghost_communication() override {}

  void solve() override {}

  void add_charge_to_field(std::size_t id, double valency) override {}
  void reset_charge_field() override {}
};

} // namespace walberla
