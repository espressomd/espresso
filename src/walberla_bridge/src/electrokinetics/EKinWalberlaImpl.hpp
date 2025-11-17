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

#include <blockforest/communication/UniformBufferedScheme.h>
#include <field/AddToStorage.h>
#include <field/FlagField.h>
#include <field/FlagUID.h>
#include <field/GhostLayerField.h>
#include <field/communication/PackInfo.h>
#include <field/vtk/FlagFieldCellFilter.h>
#include <field/vtk/VTKWriter.h>
#include <stencil/D3Q27.h>
#if defined(__CUDACC__)
#include <gpu/AddGPUFieldToStorage.h>
#include <gpu/HostFieldAllocator.h>
#include <gpu/communication/MemcpyPackInfo.h>
#include <gpu/communication/UniformGPUScheme.h>
#endif

#include "../BoundaryHandling.hpp"
#include "../BoundaryPackInfo.hpp"
#include "../utils/boundary.hpp"
#include "../utils/types_conversion.hpp"
#include "ek_kernels.hpp"
#if defined(__CUDACC__)
#include "ek_kernels.cuh"
#endif

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/BlockAndCell.hpp>
#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/EKinWalberlaBase.hpp>

#include <utils/Vector.hpp>

#include <cstddef>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

namespace walberla {

/** @brief Class that runs and controls the EK on waLBerla. */
template <std::size_t FluxCount = 13, typename FloatType = double,
          lbmpy::Arch Architecture = lbmpy::Arch::CPU>
class EKinWalberlaImpl : public EKinWalberlaBase {
  using ContinuityKernel =
      typename detail::KernelTrait<FloatType, Architecture>::ContinuityKernel;
  using DiffusiveFluxKernelUnthermalized =
      typename detail::KernelTrait<FloatType,
                                   Architecture>::DiffusiveFluxKernel;
  using DiffusiveFluxKernelThermalized = typename detail::KernelTrait<
      FloatType, Architecture>::DiffusiveFluxKernelThermalized;
  using AdvectiveFluxKernel =
      typename detail::KernelTrait<FloatType,
                                   Architecture>::AdvectiveFluxKernel;
  using FrictionCouplingKernel =
      typename detail::KernelTrait<FloatType,
                                   Architecture>::FrictionCouplingKernel;
  using DiffusiveFluxKernelElectrostaticUnthermalized =
      typename detail::KernelTrait<
          FloatType, Architecture>::DiffusiveFluxKernelElectrostatic;
  using DiffusiveFluxKernelElectrostaticThermalized =
      typename detail::KernelTrait<
          FloatType, Architecture>::DiffusiveFluxKernelElectrostaticThermalized;

  using DiffusiveFluxKernel = std::variant<DiffusiveFluxKernelUnthermalized,
                                           DiffusiveFluxKernelThermalized>;
  using DiffusiveFluxKernelElectrostatic =
      std::variant<DiffusiveFluxKernelElectrostaticUnthermalized,
                   DiffusiveFluxKernelElectrostaticThermalized>;

  using Dirichlet =
      typename detail::KernelTrait<FloatType, Architecture>::Dirichlet;
  using FixedFlux =
      typename detail::KernelTrait<FloatType, Architecture>::FixedFlux;

  using BoundaryModelDensity =
      BoundaryHandling<FloatType, FloatType, Dirichlet>;
  using BoundaryModelFlux =
      BoundaryHandling<FloatType, Vector3<FloatType>, FixedFlux>;

public:
  /** @brief Stencil for collision and streaming operations. */
  using Stencil = stencil::D3Q27;
  /** @brief Lattice model (e.g. blockforest). */
  using BlockStorage = LatticeWalberla::Lattice_T;

protected:
  template <typename FT, lbmpy::Arch AT = lbmpy::Arch::CPU> struct FieldTrait {
    // Type definitions
    using FluxField = GhostLayerField<FT, FluxCount>;
    using DensityField = GhostLayerField<FT, 1>;
    template <class Field>
    using PackInfo = field::communication::PackInfo<Field>;
    template <class Stencil>
    using RegularCommScheme =
        blockforest::communication::UniformBufferedScheme<Stencil>;
    template <class Stencil>
    using BoundaryCommScheme =
        blockforest::communication::UniformBufferedScheme<Stencil>;
  };
  using FlagField = walberla::FlagField<walberla::uint8_t>;
#if defined(__CUDACC__)
  template <typename FT> struct FieldTrait<FT, lbmpy::Arch::GPU> {
  private:
    static auto constexpr AT = lbmpy::Arch::GPU;
    template <class Field>
    using MemcpyPackInfo = gpu::communication::MemcpyPackInfo<Field>;

  public:
    template <typename Stencil>
    class UniformGPUScheme
        : public gpu::communication::UniformGPUScheme<Stencil> {
    public:
      explicit UniformGPUScheme(auto const &bf)
          : gpu::communication::UniformGPUScheme<Stencil>(
                bf, /* sendDirectlyFromGPU */ false,
                /* useLocalCommunication */ false) {}
    };
    using FluxField = gpu::GPUField<FT>;
    using DensityField = gpu::GPUField<FT>;
    template <class Field> using PackInfo = MemcpyPackInfo<Field>;
    template <class Stencil>
    using RegularCommScheme = UniformGPUScheme<Stencil>;
    template <class Stencil>
    using BoundaryCommScheme =
        blockforest::communication::UniformBufferedScheme<Stencil>;
  };
  using GPUField = gpu::GPUField<FloatType>;
#endif

  struct GhostComm {
    /** @brief Ghost communication operations. */
    enum GhostCommFlags : unsigned {
      FLB,  ///< flux boundary communication
      DENS, ///< density communication
      SIZE
    };
  };

  // "underlying" field types (`GPUField` has no f-size info at compile time)
  using _FluxField = typename FieldTrait<FloatType>::FluxField;
  using _DensityField = typename FieldTrait<FloatType>::DensityField;

public:
  using FluxField = typename FieldTrait<FloatType, Architecture>::FluxField;
  using DensityField =
      typename FieldTrait<FloatType, Architecture>::DensityField;

#if defined(__CUDACC__)
  using DensityFieldCpu = FieldTrait<FloatType, lbmpy::Arch::CPU>::DensityField;
  using FluxFieldCpu = FieldTrait<FloatType, lbmpy::Arch::CPU>::FluxField;
#endif

  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  [[nodiscard]] std::size_t stencil_size() const noexcept override {
    return FluxCount;
  }

  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same_v<FloatType, double>;
  }

private:
  FloatType m_diffusion;
  FloatType m_kT;
  FloatType m_valency;
  Utils::Vector3d m_ext_efield;
  bool m_advection;
  bool m_friction_coupling;
  unsigned int m_seed;

protected:
  // Block data access handles
  BlockDataID m_density_field_id;

  BlockDataID m_flux_field_id;

  BlockDataID m_flag_field_density_id;
  BlockDataID m_flag_field_flux_id;

#if defined(__CUDACC__)
  std::optional<BlockDataID> m_density_cpu_field_id;
  std::optional<BlockDataID> m_flux_cpu_field_id;
#endif

  /** Flag for domain cells, i.e. all cells. */
  FlagUID const Domain_flag{"domain"};
  /** Flag for boundary cells. */
  FlagUID const Boundary_flag{"boundary"};

  /** Block forest */
  std::shared_ptr<LatticeWalberla> m_lattice;

  std::unique_ptr<BoundaryModelDensity> m_boundary_density;
  std::shared_ptr<BoundaryModelFlux> m_boundary_flux;

  std::unique_ptr<DiffusiveFluxKernel> m_diffusive_flux;
  std::unique_ptr<DiffusiveFluxKernelElectrostatic>
      m_diffusive_flux_electrostatic;
  std::unique_ptr<ContinuityKernel> m_continuity;

#if defined(__CUDACC__)
  std::shared_ptr<gpu::HostFieldAllocator<FloatType>> m_host_field_allocator;
#endif

  // ResetFlux + external force
  // TODO: kernel for that
  // std::shared_ptr<ResetForce<PdfField, VectorField>> m_reset_force;

  /**
   * @brief Convenience function to add a field with a custom allocator.
   *
   * When vectorization is off, let waLBerla decide which memory allocator
   * to use. When vectorization is on, the aligned memory allocator is
   * required, otherwise <tt>cpu_vectorize_info["assume_aligned"]</tt> will
   * trigger assertions. That is because for single-precision kernels the
   * waLBerla heuristic in <tt>src/field/allocation/FieldAllocator.h</tt>
   * will fall back to @c StdFieldAlloc, yet @c AllocateAligned is needed
   * for intrinsics to work.
   */
  template <typename Field>
  auto add_to_storage(std::string const tag, FloatType value) {
    auto const &blocks = m_lattice->get_blocks();
    auto const n_ghost_layers = m_lattice->get_ghost_layers();
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      return field::addToStorage<Field>(blocks, tag, FloatType{value},
                                        field::fzyx, n_ghost_layers);
    }
#if defined(__CUDACC__)
    else {
      auto field_id = gpu::addGPUFieldToStorage<GPUField>(
          blocks, tag, Field::F_SIZE, field::fzyx, n_ghost_layers);
      if constexpr (std::is_same_v<Field, _DensityField>) {
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
          auto field = block->template getData<GPUField>(field_id);
          ek::accessor::Scalar::initialize(field, FloatType{value});
        }
      } else if constexpr (std::is_same_v<Field, _FluxField>) {
        for (auto block = blocks->begin(); block != blocks->end(); ++block) {
          auto field = block->template getData<GPUField>(field_id);
          ek::accessor::Flux::initialize(field,
                                         std::array<FloatType, FluxCount>{});
        }
      }
      return field_id;
    }
#endif
  }

  void
  reset_density_boundary_handling(std::shared_ptr<BlockStorage> const &blocks) {
    auto const [lc, uc] = m_lattice->get_local_grid_range(true);
    m_boundary_density = std::make_unique<BoundaryModelDensity>(
        blocks, m_density_field_id, m_flag_field_density_id,
        CellInterval{to_cell(lc), to_cell(uc)});
  }

  void
  reset_flux_boundary_handling(std::shared_ptr<BlockStorage> const &blocks) {
    auto const [lc, uc] = m_lattice->get_local_grid_range(true);
    m_boundary_flux = std::make_shared<BoundaryModelFlux>(
        blocks, m_flux_field_id, m_flag_field_flux_id,
        CellInterval{to_cell(lc), to_cell(uc)});
  }

  using FullCommunicator =
      typename FieldTrait<FloatType, Architecture>::template RegularCommScheme<
          typename stencil::D3Q27>;
  using BoundaryFullCommunicator =
      typename FieldTrait<FloatType, Architecture>::template BoundaryCommScheme<
          typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;
  std::shared_ptr<BoundaryFullCommunicator> m_boundary_communicator;
  std::bitset<GhostComm::SIZE> m_pending_ghost_comm;
  template <class Field>
  using PackInfo =
      typename FieldTrait<FloatType, Architecture>::template PackInfo<Field>;

public:
  EKinWalberlaImpl(std::shared_ptr<LatticeWalberla> lattice, double diffusion,
                   double kT, double valency, Utils::Vector3d const &ext_efield,
                   double density, bool advection, bool friction_coupling,
                   bool thermalized, unsigned int seed)
      : m_diffusion(FloatType_c(diffusion)), m_kT(FloatType_c(kT)),
        m_valency(FloatType_c(valency)), m_ext_efield(ext_efield),
        m_advection(advection), m_friction_coupling(friction_coupling),
        m_seed(seed), m_lattice(std::move(lattice)) {

    auto const &blocks = m_lattice->get_blocks();
    auto const n_ghost_layers = m_lattice->get_ghost_layers();

    m_density_field_id =
        add_to_storage<_DensityField>("density field", FloatType_c(density));
    m_flux_field_id =
        add_to_storage<_FluxField>("flux field", FloatType_c(0.0));

    m_continuity =
        std::make_unique<ContinuityKernel>(m_flux_field_id, m_density_field_id);

#if defined(__CUDACC__)
    m_host_field_allocator =
        std::make_shared<gpu::HostFieldAllocator<FloatType>>();
#endif

    if (thermalized) {
      set_diffusion_kernels(*m_lattice, seed);
    } else {
      set_diffusion_kernels();
    }

    // Init boundary related stuff
    m_flag_field_density_id = field::addFlagFieldToStorage<FlagField>(
        blocks, "flag field density", n_ghost_layers);
    reset_density_boundary_handling(blocks);

    m_flag_field_flux_id = field::addFlagFieldToStorage<FlagField>(
        blocks, "flag field flux", n_ghost_layers);
    reset_flux_boundary_handling(blocks);

    m_full_communication = std::make_shared<FullCommunicator>(blocks);
    m_full_communication->addPackInfo(
        std::make_shared<PackInfo<DensityField>>(m_density_field_id));
    m_boundary_communicator =
        std::make_shared<BoundaryFullCommunicator>(blocks);
    m_boundary_communicator->addPackInfo(
        std::make_shared<field::communication::BoundaryFlagPackInfo<FlagField>>(
            m_flag_field_flux_id));
    auto flux_boundary_packinfo = std::make_shared<
        field::communication::BoundaryPackInfo<FlagField, BoundaryModelFlux>>(
        m_flag_field_flux_id);
    flux_boundary_packinfo->setup_boundary_handle(m_lattice, m_boundary_flux);
    m_boundary_communicator->addPackInfo(flux_boundary_packinfo);

    m_pending_ghost_comm.set();
  }

  // Global parameters
  [[nodiscard]] double get_diffusion() const noexcept override {
    return m_diffusion;
  }
  [[nodiscard]] double get_kT() const noexcept override { return m_kT; }
  [[nodiscard]] double get_valency() const noexcept override {
    return m_valency;
  }
  [[nodiscard]] bool get_advection() const noexcept override {
    return m_advection;
  }
  [[nodiscard]] bool get_friction_coupling() const noexcept override {
    return m_friction_coupling;
  }
  [[nodiscard]] Utils::Vector3d get_ext_efield() const noexcept override {
    return m_ext_efield;
  }
  [[nodiscard]] bool is_thermalized() const noexcept override {
    return static_cast<bool>(
        std::get_if<DiffusiveFluxKernelThermalized>(&*m_diffusive_flux));
  }
  [[nodiscard]] unsigned int get_seed() const noexcept override {
    return m_seed;
  }
  [[nodiscard]] std::optional<uint64_t> get_rng_state() const override {
    auto const kernel =
        std::get_if<DiffusiveFluxKernelThermalized>(&*m_diffusive_flux);
    if (!kernel) {
      return std::nullopt;
    }
    return {static_cast<uint64_t>(kernel->getTime_step())};
  }

  void set_diffusion(double diffusion) override {
    m_diffusion = FloatType_c(diffusion);
    auto visitor = [m_diffusion = m_diffusion](auto &kernel) {
      kernel.setD(m_diffusion);
    };
    std::visit(visitor, *m_diffusive_flux);
    std::visit(visitor, *m_diffusive_flux_electrostatic);
  }

  void set_kT(double kT) override {
    m_kT = FloatType_c(kT);
    std::visit([m_kT = m_kT](auto &kernel) { kernel.setKt(m_kT); },
               *m_diffusive_flux_electrostatic);
  }

  void set_valency(double valency) override {
    m_valency = FloatType_c(valency);
    std::visit(
        [m_valency = m_valency](auto &kernel) { kernel.setZ(m_valency); },
        *m_diffusive_flux_electrostatic);
  }

  void set_advection(bool advection) override { m_advection = advection; }

  void set_friction_coupling(bool friction_coupling) override {
    m_friction_coupling = friction_coupling;
  }

  void set_rng_state(uint64_t counter) override {
    auto const kernel =
        std::get_if<DiffusiveFluxKernelThermalized>(&*m_diffusive_flux);
    auto const kernel_electrostatic =
        std::get_if<DiffusiveFluxKernelElectrostaticThermalized>(
            &*m_diffusive_flux_electrostatic);

    if (!kernel or !kernel_electrostatic) {
      throw std::runtime_error("This EK instance is unthermalized");
    }
    assert(counter <=
           static_cast<uint32_t>(std::numeric_limits<uint_t>::max()));
    kernel->setTime_step(static_cast<uint32_t>(counter));
    kernel_electrostatic->setTime_step(static_cast<uint32_t>(counter));
  }

  void set_ext_efield(Utils::Vector3d const &field) override {
    m_ext_efield = field;

    std::visit(
        [this](auto &kernel) {
          kernel.setF_ext_0(FloatType_c(m_ext_efield[0]));
          kernel.setF_ext_1(FloatType_c(m_ext_efield[1]));
          kernel.setF_ext_2(FloatType_c(m_ext_efield[2]));
        },
        *m_diffusive_flux_electrostatic);
  }

  void ghost_communication() override {
    if (m_pending_ghost_comm.test(GhostComm::DENS)) {
      (*m_full_communication)();
      m_pending_ghost_comm.reset(GhostComm::DENS);
    }
    ghost_communication_boundary();
  }

  void ghost_communication_boundary() {
    if (m_pending_ghost_comm.test(GhostComm::FLB)) {
      m_boundary_communicator->communicate();
      m_pending_ghost_comm.reset(GhostComm::FLB);
    }
  }

private:
  void set_diffusion_kernels() {
    auto kernel = DiffusiveFluxKernelUnthermalized(
        m_flux_field_id, m_density_field_id, FloatType_c(m_diffusion));
    m_diffusive_flux = std::make_unique<DiffusiveFluxKernel>(std::move(kernel));

    auto kernel_electrostatic = DiffusiveFluxKernelElectrostaticUnthermalized(
        m_flux_field_id, BlockDataID{}, m_density_field_id,
        FloatType_c(m_diffusion), FloatType_c(m_ext_efield[0]),
        FloatType_c(m_ext_efield[1]), FloatType_c(m_ext_efield[2]),
        FloatType_c(m_kT), FloatType_c(m_valency));

    m_diffusive_flux_electrostatic =
        std::make_unique<DiffusiveFluxKernelElectrostatic>(
            std::move(kernel_electrostatic));
  }

  void set_diffusion_kernels(LatticeWalberla const &lattice,
                             unsigned int seed) {
    auto const grid_dim = lattice.get_grid_dimensions();

    auto kernel = DiffusiveFluxKernelThermalized(
        m_flux_field_id, m_density_field_id, FloatType_c(m_diffusion),
        grid_dim[0], grid_dim[1], grid_dim[2], seed, 0);

    auto kernel_electrostatic = DiffusiveFluxKernelElectrostaticThermalized(
        m_flux_field_id, BlockDataID{}, m_density_field_id,
        FloatType_c(m_diffusion), FloatType_c(m_ext_efield[0]),
        FloatType_c(m_ext_efield[1]), FloatType_c(m_ext_efield[2]), grid_dim[0],
        grid_dim[1], grid_dim[2], FloatType_c(m_kT), seed, 0,
        FloatType_c(m_valency));

    auto const blocks = lattice.get_blocks();

    for (auto &block : *blocks) {
      kernel.configure(blocks, &block);
      kernel_electrostatic.configure(blocks, &block);
    }

    m_diffusive_flux = std::make_unique<DiffusiveFluxKernel>(std::move(kernel));
    m_diffusive_flux_electrostatic =
        std::make_unique<DiffusiveFluxKernelElectrostatic>(
            std::move(kernel_electrostatic));
  }

  void kernel_boundary_density() {
    for (auto &block : *m_lattice->get_blocks()) {
      (*m_boundary_density)(&block);
    }
  }

  void kernel_boundary_flux() {
    for (auto &block : *m_lattice->get_blocks()) {
      (*m_boundary_flux)(&block);
    }
  }

  void kernel_continuity() {
    for (auto &block : *m_lattice->get_blocks()) {
      (*m_continuity).run(&block);
    }
  }

  void kernel_diffusion() {
    for (auto &block : *m_lattice->get_blocks()) {
      std::visit([&block](auto &kernel) { kernel.run(&block); },
                 *m_diffusive_flux);
    }

    if (auto *kernel =
            std::get_if<DiffusiveFluxKernelThermalized>(&*m_diffusive_flux)) {
      kernel->setTime_step(kernel->getTime_step() + 1u);

      auto *kernel_electrostatic =
          std::get_if<DiffusiveFluxKernelElectrostaticThermalized>(
              &*m_diffusive_flux_electrostatic);
      kernel_electrostatic->setTime_step(kernel_electrostatic->getTime_step() +
                                         1u);
    }
  }

  void kernel_advection(std::size_t const velocity_id) {
    auto kernel = AdvectiveFluxKernel(m_flux_field_id, m_density_field_id,
                                      BlockDataID(velocity_id));
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  void kernel_friction_coupling(std::size_t const force_id,
                                double const lb_density) {
    auto kernel = FrictionCouplingKernel(
        BlockDataID(force_id), m_flux_field_id, FloatType_c(get_diffusion()),
        FloatType_c(get_kT()), FloatType(lb_density));
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  void kernel_diffusion_electrostatic(std::size_t const potential_id) {
    auto const phiID = BlockDataID(potential_id);
    std::visit([phiID](auto &kernel) { kernel.setPhiID(phiID); },
               *m_diffusive_flux_electrostatic);

    for (auto &block : *m_lattice->get_blocks()) {
      std::visit([&block](auto &kernel) { kernel.run(&block); },
                 *m_diffusive_flux_electrostatic);
    }

    if (auto *kernel_electrostatic =
            std::get_if<DiffusiveFluxKernelElectrostaticThermalized>(
                &*m_diffusive_flux_electrostatic)) {
      kernel_electrostatic->setTime_step(kernel_electrostatic->getTime_step() +
                                         1u);

      auto *kernel =
          std::get_if<DiffusiveFluxKernelThermalized>(&*m_diffusive_flux);
      kernel->setTime_step(kernel->getTime_step() + 1u);
    }
  }

  void kernel_migration() {}

  void updated_boundary_fields() {
    m_boundary_flux->boundary_update();
    m_boundary_density->boundary_update();
  }

protected:
  void integrate_vtk_writers() override {
    for (auto const &it : m_vtk_auto) {
      auto &vtk_handle = it.second;
      if (vtk_handle->enabled) {
        vtk::writeFiles(vtk_handle->ptr)();
        vtk_handle->execution_count++;
      }
    }
  }

public:
  void integrate(std::size_t potential_id, std::size_t velocity_id,
                 std::size_t force_id, double lb_density) override {

    updated_boundary_fields();

    if (get_diffusion() == 0.)
      return;

    if (get_valency() != 0.) {
      if (potential_id == walberla::BlockDataID{}) {
        throw std::runtime_error("Walberla EK: electrostatic potential enabled "
                                 "but no field accessible. potential id is " +
                                 std::to_string(potential_id));
      }
      kernel_diffusion_electrostatic(potential_id);
    } else {
      kernel_diffusion();
    }

    kernel_migration();
    kernel_boundary_flux();
    // friction coupling
    if (get_friction_coupling()) {
      if (force_id == walberla::BlockDataID{}) {
        throw std::runtime_error("Walberla EK: friction coupling enabled but "
                                 "no force field accessible. force_id is " +
                                 std::to_string(force_id) +
                                 ". Hint: LB may be inactive.");
      }
      kernel_friction_coupling(force_id, lb_density);
    }

    if (get_advection()) {
      if (velocity_id == walberla::BlockDataID{}) {
        throw std::runtime_error("Walberla EK: advection enabled but no "
                                 "velocity field accessible. velocity_id is " +
                                 std::to_string(velocity_id) +
                                 ". Hint: LB may be inactive.");
      }
      kernel_advection(velocity_id);
      kernel_boundary_flux();
    }
    kernel_continuity();

    // is this the expected behavior when reactions are included?
    kernel_boundary_density();
    m_pending_ghost_comm.set(GhostComm::DENS);

    // Handle VTK writers
    integrate_vtk_writers();
  }

  [[nodiscard]] std::size_t get_density_id() const noexcept override {
    static_assert(std::is_same_v<std::size_t, walberla::uint_t>);
    return static_cast<std::size_t>(m_density_field_id);
  }

  bool set_node_density(Utils::Vector3i const &node, double density) override {
    m_pending_ghost_comm.set(GhostComm::DENS);
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc)
      return false;

    auto density_field =
        bc->block->template getData<DensityField>(m_density_field_id);
    ek::accessor::Scalar::set(density_field, FloatType_c(density), bc->cell);
    return true;
  }

  [[nodiscard]] std::optional<double>
  get_node_density(Utils::Vector3i const &node,
                   bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);

    if (!bc)
      return std::nullopt;

    auto const density_field =
        bc->block->template getData<DensityField>(m_density_field_id);
    return {double_c(ek::accessor::Scalar::get(density_field, bc->cell))};
  }

  [[nodiscard]] std::vector<double>
  get_slice_density(Utils::Vector3i const &lower_corner,
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
          auto const density_field =
              block.template getData<DensityField>(m_density_field_id);
          auto const values = ek::accessor::Scalar::get(density_field, *bci);
          assert(values.size() == bci->numCells());
          values_size += bci->numCells();
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

  void set_slice_density(Utils::Vector3i const &lower_corner,
                         Utils::Vector3i const &upper_corner,
                         std::vector<double> const &density) override {
    m_pending_ghost_comm.set(GhostComm::DENS);
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      assert(density.size() == ci->numCells());
      for (auto &block : *lattice.get_blocks()) {
        auto const block_offset = lattice.get_block_corner(block, true);
        if (auto const bci = get_block_interval(
                lattice, lower_corner, upper_corner, block_offset, block)) {
          auto const density_field =
              block.template getData<DensityField>(m_density_field_id);
          std::vector<FloatType> values(bci->numCells());

          auto kernel = [&values, &density](unsigned const block_index,
                                            unsigned const local_index,
                                            Utils::Vector3i const &) {
            values[block_index] = numeric_cast<FloatType>(density[local_index]);
          };

          copy_block_buffer(*bci, *ci, block_offset, lower_corner, kernel);
          ek::accessor::Scalar::set(density_field, values, *bci);
        }
      }
    }
  }

  [[nodiscard]] std::optional<Utils::Vector3d>
  get_node_flux_vector(Utils::Vector3i const &node,
                       bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);

    if (!bc)
      return std::nullopt;

    auto const flux_field =
        bc->block->template getData<FluxField>(m_flux_field_id);
    return to_vector3d(ek::accessor::Flux::get_vector(flux_field, bc->cell));
  }

  std::vector<double>
  get_slice_flux_vector(Utils::Vector3i const &lower_corner,
                        Utils::Vector3i const &upper_corner) const override {
    std::vector<double> out;
    uint_t values_size = 0;
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      out = std::vector<double>(3u * ci->numCells());
      for (auto &block : *lattice.get_blocks()) {
        auto const block_offset = lattice.get_block_corner(block, true);
        if (auto const bci = get_block_interval(
                lattice, lower_corner, upper_corner, block_offset, block)) {
          auto const flux_field =
              block.template getData<FluxField>(m_flux_field_id);
          auto const values = ek::accessor::Flux::get_vector(flux_field, *bci);
          assert(values.size() == 3u * bci->numCells());
          values_size += 3u * bci->numCells();

          auto kernel = [&values, &out, this](unsigned const block_index,
                                              unsigned const local_index,
                                              Utils::Vector3i const &node) {
            if (m_boundary_flux->node_is_boundary(node)) {
              auto const &vec =
                  m_boundary_flux->get_node_value_at_boundary(node);
              for (uint_t f = 0u; f < 3u; ++f) {
                out[3u * local_index + f] = double_c(vec[f]);
              }
            } else {
              for (uint_t f = 0u; f < 3u; ++f) {
                out[3u * local_index + f] =
                    double_c(values[3u * block_index + f]);
              }
            }
          };

          copy_block_buffer(*bci, *ci, block_offset, lower_corner, kernel);
        }
      }
      assert(values_size == 3u * ci->numCells());
    }
    return out;
  }

  void clear_flux_boundaries() override {
    m_pending_ghost_comm.set(GhostComm::FLB);
    reset_flux_boundary_handling(get_lattice().get_blocks());
  }

  void clear_density_boundaries() override {
    reset_density_boundary_handling(get_lattice().get_blocks());
  }

  bool set_node_flux_boundary(Utils::Vector3i const &node,
                              Utils::Vector3d const &flux) override {
    m_pending_ghost_comm.set(GhostComm::FLB);
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_flux->set_node_value_at_boundary(
        node, to_vector3<FloatType>(flux), *bc);
    return true;
  }

  [[nodiscard]] std::optional<Utils::Vector3d>
  get_node_flux_at_boundary(Utils::Vector3i const &node,
                            bool consider_ghosts = false) const override {
    assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::FLB)));
    auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc or !m_boundary_flux->node_is_boundary(node))
      return std::nullopt;

    return {to_vector3d(m_boundary_flux->get_node_value_at_boundary(node))};
  }

  bool remove_node_from_flux_boundary(Utils::Vector3i const &node) override {
    m_pending_ghost_comm.set(GhostComm::FLB);
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_flux->remove_node_from_boundary(node, *bc);
    return true;
  }

  bool set_node_density_boundary(Utils::Vector3i const &node,
                                 double density) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_density->set_node_value_at_boundary(node, FloatType_c(density),
                                                   *bc);

    return true;
  }

  [[nodiscard]] std::optional<double>
  get_node_density_at_boundary(Utils::Vector3i const &node,
                               bool consider_ghosts = false) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc or !m_boundary_density->node_is_boundary(node))
      return std::nullopt;

    return {double_c(m_boundary_density->get_node_value_at_boundary(node))};
  }

  void set_slice_density_boundary(
      Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
      std::vector<std::optional<double>> const &density) override {
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      auto const local_offset = std::get<0>(lattice.get_local_grid_range());
      auto const lower_cell = ci->min();
      auto const upper_cell = ci->max();
      auto it = density.begin();
      assert(density.size() == ci->numCells());
      for (auto x = lower_cell.x(); x <= upper_cell.x(); ++x) {
        for (auto y = lower_cell.y(); y <= upper_cell.y(); ++y) {
          for (auto z = lower_cell.z(); z <= upper_cell.z(); ++z) {
            auto const node = local_offset + Utils::Vector3i{{x, y, z}};
            auto const bc = get_block_and_cell(lattice, node, false);
            auto const &opt = *it;
            if (opt) {
              m_boundary_density->set_node_value_at_boundary(
                  node, FloatType_c(*opt), *bc);
            } else {
              m_boundary_density->remove_node_from_boundary(node, *bc);
            }
            ++it;
          }
        }
      }
    }
  }

  [[nodiscard]] std::vector<std::optional<double>>
  get_slice_density_at_boundary(
      Utils::Vector3i const &lower_corner,
      Utils::Vector3i const &upper_corner) const override {
    std::vector<std::optional<double>> out;
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      auto const local_offset = std::get<0>(lattice.get_local_grid_range());
      auto const lower_cell = ci->min();
      auto const upper_cell = ci->max();
      auto const n_values = ci->numCells();
      out.reserve(n_values);
      for (auto x = lower_cell.x(); x <= upper_cell.x(); ++x) {
        for (auto y = lower_cell.y(); y <= upper_cell.y(); ++y) {
          for (auto z = lower_cell.z(); z <= upper_cell.z(); ++z) {
            auto const node = local_offset + Utils::Vector3i{{x, y, z}};
            if (m_boundary_density->node_is_boundary(node)) {
              out.emplace_back(double_c(
                  m_boundary_density->get_node_value_at_boundary(node)));
            } else {
              out.emplace_back(std::nullopt);
            }
          }
        }
      }
      assert(out.size() == n_values);
    }
    return out;
  }

  void set_slice_flux_boundary(
      Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
      std::vector<std::optional<Utils::Vector3d>> const &flux) override {
    m_pending_ghost_comm.set(GhostComm::FLB);
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      auto const local_offset = std::get<0>(lattice.get_local_grid_range());
      auto const lower_cell = ci->min();
      auto const upper_cell = ci->max();
      auto it = flux.begin();
      assert(flux.size() == ci->numCells());
      for (auto x = lower_cell.x(); x <= upper_cell.x(); ++x) {
        for (auto y = lower_cell.y(); y <= upper_cell.y(); ++y) {
          for (auto z = lower_cell.z(); z <= upper_cell.z(); ++z) {
            auto const node = local_offset + Utils::Vector3i{{x, y, z}};
            auto const bc = get_block_and_cell(lattice, node, false);
            auto const &opt = *it;
            if (opt) {
              m_boundary_flux->set_node_value_at_boundary(
                  node, to_vector3<FloatType>(*opt), *bc);
            } else {
              m_boundary_flux->remove_node_from_boundary(node, *bc);
            }
            ++it;
          }
        }
      }
    }
  }

  [[nodiscard]] std::vector<std::optional<Utils::Vector3d>>
  get_slice_flux_at_boundary(
      Utils::Vector3i const &lower_corner,
      Utils::Vector3i const &upper_corner) const override {
    std::vector<std::optional<Utils::Vector3d>> out;
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      auto const local_offset = std::get<0>(lattice.get_local_grid_range());
      auto const lower_cell = ci->min();
      auto const upper_cell = ci->max();
      auto const n_values = ci->numCells();
      out.reserve(n_values);
      for (auto x = lower_cell.x(); x <= upper_cell.x(); ++x) {
        for (auto y = lower_cell.y(); y <= upper_cell.y(); ++y) {
          for (auto z = lower_cell.z(); z <= upper_cell.z(); ++z) {
            auto const node = local_offset + Utils::Vector3i{{x, y, z}};
            if (m_boundary_flux->node_is_boundary(node)) {
              out.emplace_back(to_vector3d(
                  m_boundary_flux->get_node_value_at_boundary(node)));
            } else {
              out.emplace_back(std::nullopt);
            }
          }
        }
      }
      assert(out.size() == n_values);
    }
    return out;
  }

  [[nodiscard]] std::vector<bool>
  get_slice_is_boundary(Utils::Vector3i const &lower_corner,
                        Utils::Vector3i const &upper_corner) const override {
    std::vector<bool> out;
    auto const &lattice = get_lattice();
    if (auto const ci = get_interval(lattice, lower_corner, upper_corner)) {
      auto const local_offset = std::get<0>(lattice.get_local_grid_range());
      auto const lower_cell = ci->min();
      auto const upper_cell = ci->max();
      auto const n_values = ci->numCells();
      out.reserve(n_values);
      for (auto x = lower_cell.x(); x <= upper_cell.x(); ++x) {
        for (auto y = lower_cell.y(); y <= upper_cell.y(); ++y) {
          for (auto z = lower_cell.z(); z <= upper_cell.z(); ++z) {
            auto const node = local_offset + Utils::Vector3i{{x, y, z}};
            out.emplace_back(m_boundary_density->node_is_boundary(node) or
                             m_boundary_flux->node_is_boundary(node));
          }
        }
      }
      assert(out.size() == n_values);
    }
    return out;
  }

  bool remove_node_from_density_boundary(Utils::Vector3i const &node) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_density->remove_node_from_boundary(node, *bc);

    return true;
  }

  [[nodiscard]] std::optional<bool>
  get_node_is_flux_boundary(Utils::Vector3i const &node,
                            bool consider_ghosts) const override {
    assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::FLB)));
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return std::nullopt;

    return {m_boundary_flux->node_is_boundary(node)};
  }

  [[nodiscard]] std::optional<bool>
  get_node_is_density_boundary(Utils::Vector3i const &node,
                               bool consider_ghosts) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return std::nullopt;

    return {m_boundary_density->node_is_boundary(node)};
  }

  [[nodiscard]] std::optional<bool>
  get_node_is_boundary(Utils::Vector3i const &node,
                       bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return std::nullopt;

    return {m_boundary_density->node_is_boundary(node) or
            m_boundary_flux->node_is_boundary(node)};
  }

  void update_flux_boundary_from_shape(
      const std::vector<int> &raster_flat,
      const std::vector<double> &data_flat) override {
    m_pending_ghost_comm.set(GhostComm::FLB);
    auto const grid_size = get_lattice().get_grid_dimensions();
    auto const data = fill_3D_vector_array(data_flat, grid_size);
    set_boundary_from_grid(*m_boundary_flux, get_lattice(), raster_flat, data);
    reallocate_flux_boundary_field();
  }

  void update_density_boundary_from_shape(
      const std::vector<int> &raster_flat,
      const std::vector<double> &data_flat) override {
    auto const grid_size = get_lattice().get_grid_dimensions();
    auto const data = fill_3D_scalar_array(data_flat, grid_size);
    set_boundary_from_grid(*m_boundary_density, get_lattice(), raster_flat,
                           data);
    reallocate_density_boundary_field();
  }

  void reallocate_flux_boundary_field() { m_boundary_flux->boundary_update(); }

  void reallocate_density_boundary_field() {
    m_boundary_density->boundary_update();
  }

  [[nodiscard]] LatticeWalberla const &get_lattice() const noexcept override {
    return *m_lattice;
  }

  [[nodiscard]] bool is_gpu() const noexcept override {
    return Architecture == lbmpy::Arch::GPU;
  }

  void register_vtk_field_filters(walberla::vtk::VTKOutput &vtk_obj) override {
    field::FlagFieldCellFilter<FlagField> fluid_filter(m_flag_field_density_id);
    fluid_filter.addFlag(Boundary_flag);
    vtk_obj.addCellExclusionFilter(fluid_filter);
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

#if defined(__CUDACC__)
  template <typename OutputType = float>
  class DensityVTKWriter : public VTKWriter<DensityFieldCpu, 1u, OutputType> {
  public:
    using Base = VTKWriter<DensityFieldCpu, 1u, OutputType>;
    using Base::Base;
    using Base::evaluate;

  protected:
    OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                        cell_idx_t const z, cell_idx_t const) override {
      WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
      auto const density = ek::accessor::Scalar::get(this->m_field, {x, y, z});
      return numeric_cast<OutputType>(this->m_conversion * density);
    }
  };
#else
  template <typename OutputType = float>
  class DensityVTKWriter : public VTKWriter<DensityField, 1u, OutputType> {
  public:
    using Base = VTKWriter<DensityField, 1u, OutputType>;
    using Base::Base;
    using Base::evaluate;

  protected:
    OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                        cell_idx_t const z, cell_idx_t const) override {
      WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
      auto const density = ek::accessor::Scalar::get(this->m_field, {x, y, z});
      return numeric_cast<OutputType>(this->m_conversion * density);
    }
  };
#endif

#if defined(__CUDACC__)
  template <typename OutputType = float>
  class FluxVTKWriter : public VTKWriter<FluxFieldCpu, 3u, OutputType> {
  public:
    using Base = VTKWriter<FluxFieldCpu, 3u, OutputType>;
    using Base::Base;
    using Base::evaluate;

  protected:
    OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                        cell_idx_t const z, cell_idx_t const f) override {
      WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
      auto const flux =
          ek::accessor::Flux::get_vector(this->m_field, {x, y, z});
      return numeric_cast<OutputType>(this->m_conversion * flux[uint_c(f)]);
    }
  };
#else
  template <typename OutputType = float>
  class FluxVTKWriter : public VTKWriter<FluxField, 3u, OutputType> {
  public:
    using Base = VTKWriter<FluxField, 3u, OutputType>;
    using Base::Base;
    using Base::evaluate;

  protected:
    OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                        cell_idx_t const z, cell_idx_t const f) override {
      WALBERLA_ASSERT_NOT_NULLPTR(this->m_field);
      auto const flux =
          ek::accessor::Flux::get_vector(this->m_field, {x, y, z});
      return numeric_cast<OutputType>(this->m_conversion * flux[uint_c(f)]);
    }
  };
#endif

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
                m_lattice->get_ghost_layers(), m_host_field_allocator);
          }
        };
#endif
    if (flag_observables & static_cast<int>(EKOutputVTK::density)) {
      auto const unit_conversion = FloatType_c(units.at("density"));
#if defined(__CUDACC__)
      if constexpr (Architecture == lbmpy::Arch::GPU) {
        auto const &blocks = m_lattice->get_blocks();
        allocate_cpu_field_if_empty.template operator()<DensityFieldCpu>(
            blocks, "density_cpu", m_density_cpu_field_id);
        vtk_obj.addBeforeFunction(
            gpu::fieldCpyFunctor<DensityFieldCpu, DensityField>(
                blocks, *m_density_cpu_field_id, m_density_field_id));
        vtk_obj.addCellDataWriter(make_shared<DensityVTKWriter<float>>(
            *m_density_cpu_field_id, "density", unit_conversion));
      } else {
#endif
        vtk_obj.addCellDataWriter(make_shared<DensityVTKWriter<float>>(
            m_density_field_id, "density", unit_conversion));
#if defined(__CUDACC__)
      }
#endif
    }
    if (flag_observables & static_cast<int>(EKOutputVTK::flux)) {
      auto const unit_conversion = FloatType_c(units.at("flux"));
#if defined(__CUDACC__)
      if constexpr (Architecture == lbmpy::Arch::GPU) {
        auto const &blocks = m_lattice->get_blocks();
        allocate_cpu_field_if_empty.template operator()<FluxFieldCpu>(
            blocks, "flux_cpu", m_flux_cpu_field_id);
        vtk_obj.addBeforeFunction(gpu::fieldCpyFunctor<FluxFieldCpu, FluxField>(
            blocks, *m_flux_cpu_field_id, m_flux_field_id));
        vtk_obj.addCellDataWriter(make_shared<FluxVTKWriter<float>>(
            *m_flux_cpu_field_id, "flux", unit_conversion));
      } else {
#endif
        vtk_obj.addCellDataWriter(make_shared<FluxVTKWriter<float>>(
            m_flux_field_id, "flux", unit_conversion));
#if defined(__CUDACC__)
      }
#endif
    }
  }

  ~EKinWalberlaImpl() override = default;
};

} // namespace walberla
