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
#include <field/communication/PackInfo.h>
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
#include <gpu/communication/MemcpyPackInfo.h>
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

#include <heffte.h>
#include <heffte_backends.h>
#include <heffte_geometry.h>

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

template <typename T, std::size_t N>
auto to_array(Utils::Vector<T, N> const &vec) {
  std::array<T, N> res{};
  std::ranges::copy(vec, res.begin());
  return res;
}

inline int pos_to_linear_index(int x, int y, int z, auto const &dim) {
  return (z * dim[1] + y) * dim[0] + x;
}

template <typename FloatType, lbmpy::Arch Architecture>
class PoissonSolverFFT : public PoissonSolver {
private:
  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

protected:
  template <typename FT, lbmpy::Arch AT = lbmpy::Arch::CPU> struct FieldTrait {
    using ComplexType = std::complex<FloatType>;
    using PotentialField = field::GhostLayerField<FT, 1u>;
    using PotentialFieldVTK = field::GhostLayerField<FT, 1u>;
    template <class Field>
    using PackInfo = field::communication::PackInfo<Field>;
    template <class Stencil>
    using FullCommunicator =
        blockforest::communication::UniformBufferedScheme<Stencil>;
  };

#if defined(__CUDACC__)
  template <typename FT> struct FieldTrait<FT, lbmpy::Arch::GPU> {
    using ComplexType = std::conditional_t<std::is_same_v<FloatType, float>,
                                           cufftComplex, cufftDoubleComplex>;
    using PotentialField = gpu::GPUField<FT>;
    using PotentialFieldVTK = field::GhostLayerField<FT, 1u>;
    template <class Field>
    using PackInfo = gpu::communication::MemcpyPackInfo<Field>;
    template <class Stencil>
    using FullCommunicator = gpu::communication::UniformGPUScheme<Stencil>;
  };
#endif

public:
  using ComplexType = FieldTrait<FloatType, Architecture>::ComplexType;
  using PotentialField = FieldTrait<FloatType, Architecture>::PotentialField;
  using PotentialFieldVTK =
      FieldTrait<FloatType, Architecture>::PotentialFieldVTK;
  using FullCommunicator =
      FieldTrait<FloatType,
                 Architecture>::template FullCommunicator<stencil::D3Q27>;
  template <class Field>
  using PackInfo =
      FieldTrait<FloatType, Architecture>::template PackInfo<Field>;

protected:
  template <typename ComplexType> struct heffte_container {
#if defined(__CUDACC__)
    using backend = heffte::backend::cufft;
#else  // __CUDACC__
    using backend = heffte::backend::fftw;
#endif // __CUDACC__
    std::unique_ptr<heffte::box3d<>> box_in;
    std::unique_ptr<heffte::box3d<>> box_out;
    std::unique_ptr<heffte::fft3d<backend>> fft;
    std::unique_ptr<heffte::fft3d<backend>::buffer_container<ComplexType>>
        buffer;

    heffte_container(heffte::plan_options const &options,
                     auto const &grid_range) {
      auto const offset_vec = Utils::Vector3i({1, 1, 1});
      auto const order = Utils::Vector3i({0, 1, 2});
      box_in = std::make_unique<heffte::box3d<>>(
          to_array(Utils::Vector3i(grid_range.first)),
          to_array(grid_range.second - offset_vec), to_array(order));
      box_out = std::make_unique<heffte::box3d<>>(
          to_array(Utils::Vector3i(grid_range.first)),
          to_array(grid_range.second - offset_vec), to_array(order));
      fft = std::make_unique<heffte::fft3d<backend>>(*box_in, *box_out,
                                                     MPI_COMM_WORLD, options);
      buffer = std::make_unique<
          heffte::fft3d<backend>::buffer_container<ComplexType>>(
          fft->size_workspace());
    }
  };

private:
  BlockDataID m_potential_field_with_ghosts_id;
#if defined(__CUDACC__)
  BlockDataID m_potential_field_id;
  BlockDataID m_greens_function_field_id;
  BlockDataID m_potential_fourier_id;
#endif // __CUDACC__

  std::unique_ptr<heffte_container<ComplexType>> heffte;
  std::shared_ptr<FullCommunicator> m_full_communication;

#if defined(__CUDACC__)
  std::optional<BlockDataID> m_potential_cpu_field_id;
  std::shared_ptr<gpu::HostFieldAllocator<FloatType>> m_host_field_allocator;

  using GreenFunctionField = gpu::GPUField<FloatType>;
  using PotentialFourier = gpu::GPUField<ComplexType>;

  walberla::gpu::Kernel<void (*)(walberla::gpu::FieldAccessor<ComplexType>,
                                 walberla::gpu::FieldAccessor<FloatType>)>
      kernel_greens;
  walberla::gpu::Kernel<void (*)(walberla::gpu::FieldAccessor<FloatType>,
                                 walberla::gpu::FieldAccessor<FloatType>)>
      kernel_move_fields;
#else  // __CUDACC__
  std::vector<FloatType> m_greens;
  std::vector<FloatType> m_potential;
  std::vector<ComplexType> m_potential_fourier;
#endif // __CUDACC__

public:
  ~PoissonSolverFFT() override = default;
  PoissonSolverFFT(std::shared_ptr<LatticeWalberla> lattice,
                   double permittivity)
      : PoissonSolver(std::move(lattice), permittivity)
#if defined(__CUDACC__)
        ,
        kernel_greens(gpu::make_kernel(
            multiply_by_greens_function<FloatType, ComplexType>)),
        kernel_move_fields(gpu::make_kernel(move_field<FloatType>))
#endif
  {
    auto blocks = get_lattice().get_blocks();
#if defined(__CUDACC__)
    m_potential_field_id = gpu::addGPUFieldToStorage<PotentialField>(
        blocks, "potential field", 1u, field::fzyx, 0u, false);
    m_potential_field_with_ghosts_id =
        gpu::addGPUFieldToStorage<PotentialField>(
            blocks, "potential field with ghosts", 1u, field::fzyx,
            get_lattice().get_ghost_layers());
    m_greens_function_field_id = gpu::addGPUFieldToStorage<GreenFunctionField>(
        blocks, "greens function", 1u, field::fzyx, 0u, false);
    m_potential_fourier_id = gpu::addGPUFieldToStorage<PotentialFourier>(
        blocks, "fourier field", 1u, field::fzyx, 0u, false);

    auto const grid_range = get_lattice().get_local_grid_range();
    auto const global_dim = get_lattice().get_grid_dimensions();
    for (auto &block : *blocks) {
      auto green_field = block.template getData<GreenFunctionField>(
          m_greens_function_field_id);
      auto kernel = gpu::make_kernel(create_greens_function<FloatType>);
      kernel.addFieldIndexingParam(
          gpu::FieldIndexing<FloatType>::xyz(*green_field));
      kernel.addParam(grid_range.first[0]);
      kernel.addParam(grid_range.first[1]);
      kernel.addParam(grid_range.first[2]);
      kernel.addParam(grid_range.second[0]);
      kernel.addParam(grid_range.second[1]);
      kernel.addParam(grid_range.second[2]);
      kernel.addParam(global_dim[0]);
      kernel.addParam(global_dim[1]);
      kernel.addParam(global_dim[2]);
      kernel();

      auto potential =
          block.template getData<PotentialField>(m_potential_field_id);
      auto potential_ghosts = block.template getData<PotentialField>(
          m_potential_field_with_ghosts_id);
      auto green = block.template getData<GreenFunctionField>(
          m_greens_function_field_id);
      auto fourier =
          block.template getData<PotentialFourier>(m_potential_fourier_id);

      kernel_greens =
          gpu::make_kernel(multiply_by_greens_function<FloatType, ComplexType>);
      kernel_greens.addFieldIndexingParam(
          gpu::FieldIndexing<ComplexType>::allInner(*fourier));
      kernel_greens.addFieldIndexingParam(
          gpu::FieldIndexing<FloatType>::allInner(*green));

      kernel_move_fields = gpu::make_kernel(move_field<FloatType>);
      kernel_move_fields.addFieldIndexingParam(
          gpu::FieldIndexing<FloatType>::xyz(*potential_ghosts));
      kernel_move_fields.addFieldIndexingParam(
          gpu::FieldIndexing<FloatType>::xyz(*potential));
    }

    m_host_field_allocator =
        std::make_shared<gpu::HostFieldAllocator<FloatType>>();
#else  // __CUDACC__
    m_potential_field_with_ghosts_id = field::addToStorage<PotentialField>(
        blocks, "potential field with ghosts", 0., field::fzyx,
        get_lattice().get_ghost_layers());
#endif // __CUDACC__

    m_full_communication = std::make_shared<FullCommunicator>(blocks);
    m_full_communication->addPackInfo(
        std::make_shared<PackInfo<PotentialField>>(
            m_potential_field_with_ghosts_id));
  }

  [[nodiscard]] bool is_gpu() const noexcept override {
    return Architecture == lbmpy::Arch::GPU;
  }

  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same_v<FloatType, double>;
  }

  std::size_t get_potential_field_id() const noexcept override {
    return static_cast<std::size_t>(m_potential_field_with_ghosts_id);
  }

  void setup_fft([[maybe_unused]] bool use_gpu_aware) override {
    auto const grid_range = get_lattice().get_local_grid_range();
    heffte::plan_options options = heffte::default_options<
        typename heffte_container<ComplexType>::backend>();
#if defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      options.use_reorder = false;
      options.algorithm = heffte::reshape_algorithm::p2p_plined;
      options.use_pencils = true;
      options.use_gpu_aware = use_gpu_aware;
    }
#endif
    heffte =
        std::make_unique<heffte_container<ComplexType>>(options, grid_range);
#if not defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      m_potential = std::vector<FloatType>(heffte->fft->size_inbox());
      m_greens = std::vector<FloatType>(heffte->fft->size_outbox());
      m_potential_fourier =
          std::vector<ComplexType>(heffte->fft->size_outbox());
      auto const dim = grid_range.second - grid_range.first;
      auto const global_dim = get_lattice().get_grid_dimensions();
      for (int x = 0; x < dim[0]; x++) {
        for (int y = 0; y < dim[1]; y++) {
          for (int z = 0; z < dim[2]; z++) {
            m_greens[pos_to_linear_index(x, y, z, dim)] =
                greens_function<FloatType>(x + grid_range.first[0],
                                           y + grid_range.first[1],
                                           z + grid_range.first[2], global_dim);
          }
        }
      }
    }
#endif
    reset_charge_field();
  }

  [[nodiscard]] std::optional<double>
  get_node_potential(Utils::Vector3i const &node,
                     bool consider_ghosts = false) override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);

    if (not bc or get_potential_field_id() == 0u)
      return std::nullopt;

    auto const potential_field = bc->block->template getData<PotentialField>(
        m_potential_field_with_ghosts_id);
    return {double_c(
        walberla::ek::accessor::Scalar::get(potential_field, bc->cell))};
  }

  [[nodiscard]] std::vector<double>
  get_slice_potential(Utils::Vector3i const &lower_corner,
                      Utils::Vector3i const &upper_corner) const override {
    std::vector<double> out;
    uint_t values_size{0u};
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      out = std::vector<double>(ci->numCells());
      for (auto &block : *lattice.get_blocks()) {
        auto const block_offset = lattice.get_block_corner(block, true);
        if (auto const bci = get_block_interval(
                lattice, lower_corner, upper_corner, block_offset, block)) {
          auto const potential_field = block.template getData<PotentialField>(
              m_potential_field_with_ghosts_id);
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

  void solve() override {
#if not defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      auto grid_range = get_lattice().get_local_grid_range();
      auto dim = grid_range.second - grid_range.first;
      heffte->fft->forward(m_potential.data(), m_potential_fourier.data(),
                           heffte->buffer->data());
      std::ranges::transform(m_potential_fourier, m_greens,
                             m_potential_fourier.begin(), std::multiplies<>{});
      heffte->fft->backward(m_potential_fourier.data(), m_potential.data(),
                            heffte->buffer->data());

      for (auto &block : *get_lattice().get_blocks()) {
        auto potential_with_ghosts = block.template getData<PotentialField>(
            m_potential_field_with_ghosts_id);
        for (int x = 0; x < dim[0]; x++) {
          for (int y = 0; y < dim[1]; y++) {
            for (int z = 0; z < dim[2]; z++) {
              potential_with_ghosts->get(x, y, z) =
                  m_potential[pos_to_linear_index(x, y, z, dim)];
            }
          }
        }
      }
    }
#endif
#if defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      for (auto &block : *get_lattice().get_blocks()) {
        auto potential =
            block.template getData<PotentialField>(m_potential_field_id);
        auto fourier =
            block.template getData<PotentialFourier>(m_potential_fourier_id);
        FloatType *_data_potential = potential->dataAt(0, 0, 0, 0);
        ComplexType *_data_fourier = fourier->dataAt(0, 0, 0, 0);
        heffte->fft->forward(_data_potential, _data_fourier,
                             heffte->buffer->data());
        kernel_greens();
        heffte->fft->backward(_data_fourier, _data_potential,
                              heffte->buffer->data());
        kernel_move_fields();
      }
    }
#endif
    ghost_communication();
    integrate_vtk_writers();
  }

  void add_charge_to_field(std::size_t id, double valency,
                           bool is_double_precision) override {
    auto const factor = FloatType_c(valency / get_permittivity());
    auto const density_id = BlockDataID(id);
#if not defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      auto grid_range = get_lattice().get_local_grid_range();
      auto dim = grid_range.second - grid_range.first;
      for (auto &block : *get_lattice().get_blocks()) {
        auto density_field = block.template getData<PotentialField>(density_id);
        for (int x = 0; x < dim[0]; x++) {
          for (int y = 0; y < dim[1]; y++) {
            for (int z = 0; z < dim[2]; z++) {
              m_potential[pos_to_linear_index(x, y, z, dim)] +=
                  factor * density_field->get(x, y, z);
            }
          }
        }
      }
    }
#endif
#if defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      for (auto &block : *get_lattice().get_blocks()) {
        auto field =
            block.template getData<PotentialField>(m_potential_field_id);
        auto density_field =
            block.template getData<gpu::GPUField<FloatType>>(density_id);
        add_fields(field, density_field, FloatType_c(factor));
      }
    }
#endif
  }

  void reset_charge_field() override {
#if not defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      auto const grid_range = get_lattice().get_local_grid_range();
      auto const dim = grid_range.second - grid_range.first;
      for (int x = 0; x < dim[0]; x++) {
        for (int i = 0; i < m_potential_fourier.size(); i++) {
          m_potential_fourier[i] *= m_greens[i];
        }
        for (int y = 0; y < dim[1]; y++) {
          for (int z = 0; z < dim[2]; z++) {
            m_potential[pos_to_linear_index(x, y, z, dim)] = FloatType(0.0);
          }
        }
      }
    }
#endif
#if defined(__CUDACC__)
    if constexpr (Architecture == lbmpy::Arch::GPU) {
      // the FFT-solver re-uses the potential field for the charge
      for (auto &block : *get_lattice().get_blocks()) {
        auto field =
            block.template getData<PotentialField>(m_potential_field_id);
        ek::accessor::Scalar::initialize(field, FloatType_c(0.));
      }
    }
#endif
  }

protected:
  void ghost_communication() { (*m_full_communication)(); }

  void integrate_vtk_writers() override {
    for (auto const &it : m_vtk_auto) {
      auto &vtk_handle = it.second;
      if (vtk_handle->enabled) {
        vtk::writeFiles(vtk_handle->ptr)();
        vtk_handle->execution_count++;
      }
    }
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
      : public VTKWriter<PotentialFieldVTK, 1u, OutputType> {
  public:
    using Base = VTKWriter<PotentialFieldVTK, 1u, OutputType>;
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

public:
  void register_vtk_field_writers(walberla::vtk::VTKOutput &vtk_obj,
                                  LatticeModel::units_map const &units,
                                  int flag_observables) override {
#if defined(__CUDACC__)
    auto const allocate_cpu_field_if_empty =
        [&]<typename Field>(auto const &blocks, std::string name,
                            std::optional<BlockDataID> &cpu_field) {
          if (not cpu_field) {
            cpu_field = field::addToStorage<Field>(
                blocks, name, FloatType{0}, field::fzyx,
                get_lattice().get_ghost_layers(), m_host_field_allocator);
          }
        };
#endif
    if (flag_observables & static_cast<int>(EKPoissonOutputVTK::potential)) {
      auto const unit_conversion = FloatType_c(units.at("potential"));
#if defined(__CUDACC__)
      if constexpr (Architecture == lbmpy::Arch::GPU) {
        auto const &blocks = get_lattice().get_blocks();
        allocate_cpu_field_if_empty.template operator()<PotentialFieldVTK>(
            blocks, "potential_cpu", m_potential_cpu_field_id);
        vtk_obj.addBeforeFunction(
            gpu::fieldCpyFunctor<PotentialFieldVTK, PotentialField>(
                blocks, *m_potential_cpu_field_id,
                m_potential_field_with_ghosts_id));
        vtk_obj.addCellDataWriter(make_shared<PotentialVTKWriter<float>>(
            *m_potential_cpu_field_id, "potential", unit_conversion));
      }
#endif
      if constexpr (Architecture == lbmpy::Arch::CPU) {
        vtk_obj.addCellDataWriter(make_shared<PotentialVTKWriter<float>>(
            m_potential_field_with_ghosts_id, "potential", unit_conversion));
      }
    }
  }

private:
#if defined(__CUDACC__)
  void add_fields(PotentialField *field_out,
                  gpu::GPUField<FloatType> *field_add, FloatType factor) {
    auto kernel = gpu::make_kernel(add_fields_with_factor<FloatType>);
    kernel.addFieldIndexingParam(
        gpu::FieldIndexing<FloatType>::xyz(*field_out));
    kernel.addFieldIndexingParam(
        gpu::FieldIndexing<FloatType>::xyz(*field_add));
    kernel.addParam(factor);
    kernel();
  }
#endif
};

} // namespace walberla
