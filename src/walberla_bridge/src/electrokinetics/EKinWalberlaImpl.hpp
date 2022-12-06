/*
 * Copyright (C) 2022 The ESPResSo project
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
#include <field/GhostLayerField.h>
#include <field/communication/PackInfo.h>
#include <field/vtk/FlagFieldCellFilter.h>
#include <field/vtk/VTKWriter.h>
#include <lbm/lattice_model/D3Q27.h>
#include <timeloop/SweepTimeloop.h>

#include "walberla_bridge/BlockAndCell.hpp"
#include "walberla_bridge/LatticeWalberla.hpp"
#include "walberla_bridge/electrokinetics/EKinWalberlaBase.hpp"
#include "walberla_bridge/utils/boundary_utils.hpp"
#include "walberla_bridge/utils/walberla_utils.hpp"

#include "ek_kernels.hpp"

#include <utils/Vector.hpp>

#include <cstddef>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "../BoundaryHandling.hpp"

namespace walberla {

/** @brief Class that runs and controls the EK on waLBerla. */
template <size_t FluxCount = 13, typename FloatType = double>
class EKinWalberlaImpl : public EKinWalberlaBase {
  template <typename T> inline FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  using ContinuityKernel =
      typename detail::KernelTrait<FloatType>::ContinuityKernel;
  using DiffusiveFluxKernel =
      typename detail::KernelTrait<FloatType>::DiffusiveFluxKernel;
  using AdvectiveFluxKernel =
      typename detail::KernelTrait<FloatType>::AdvectiveFluxKernel;
  using FrictionCouplingKernel =
      typename detail::KernelTrait<FloatType>::FrictionCouplingKernel;
  using DiffusiveFluxKernelElectrostatic =
      typename detail::KernelTrait<FloatType>::DiffusiveFluxKernelElectrostatic;

  using Dirichlet = typename detail::KernelTrait<FloatType>::Dirichlet;
  using FixedFlux = typename detail::KernelTrait<FloatType>::FixedFlux;

protected:
  // Type definitions
  using FluxField = GhostLayerField<FloatType, FluxCount>;
  using FlagField = walberla::FlagField<walberla::uint8_t>;
  using DensityField = GhostLayerField<FloatType, 1>;

  using BoundaryModelDensity = BoundaryHandling<FloatType, Dirichlet>;
  using BoundaryModelFlux = BoundaryHandling<Vector3<FloatType>, FixedFlux>;

private:
  FloatType m_diffusion;
  FloatType m_kT;
  FloatType m_valency;
  Utils::Vector3d m_ext_efield;
  bool m_advection;
  bool m_friction_coupling;

protected:
  // Block data access handles
  BlockDataID m_density_field_id;
  BlockDataID m_density_field_flattened_id;

  BlockDataID m_flux_field_id;
  BlockDataID m_flux_field_flattened_id;

  BlockDataID m_flag_field_density_id;
  BlockDataID m_flag_field_flux_id;

  /** Block forest */
  std::shared_ptr<LatticeWalberla> m_lattice;

  std::unique_ptr<BoundaryModelDensity> m_boundary_density;
  std::unique_ptr<BoundaryModelFlux> m_boundary_flux;

  std::unique_ptr<ContinuityKernel> m_continuity;

  // ResetFlux + external force
  // TODO: kernel for that
  // std::shared_ptr<ResetForce<PdfField, VectorField>> m_reset_force;

  [[nodiscard]] std::size_t stencil_size() const override { return FluxCount; }

  void reset_density_boundary_handling() {
    auto const &blocks = get_lattice().get_blocks();
    m_boundary_density = std::make_unique<BoundaryModelDensity>(
        blocks, m_density_field_id, m_flag_field_density_id);
  }

  void reset_flux_boundary_handling() {
    auto const &blocks = get_lattice().get_blocks();
    m_boundary_flux = std::make_unique<BoundaryModelFlux>(
        blocks, m_flux_field_id, m_flag_field_flux_id);
  }

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;

public:
  EKinWalberlaImpl(std::shared_ptr<LatticeWalberla> lattice, double diffusion,
                   double kT, double valency, const Utils::Vector3d &ext_efield,
                   double density, bool advection, bool friction_coupling)
      : m_diffusion(FloatType_c(diffusion)), m_kT(FloatType_c(kT)),
        m_valency(FloatType_c(valency)), m_ext_efield(ext_efield),
        m_advection(advection), m_friction_coupling(friction_coupling),
        m_lattice(std::move(lattice)) {
    m_density_field_id = field::addToStorage<DensityField>(
        m_lattice->get_blocks(), "density field", FloatType_c(density),
        field::fzyx, m_lattice->get_ghost_layers());
    m_density_field_flattened_id =
        field::addFlattenedShallowCopyToStorage<DensityField>(
            m_lattice->get_blocks(), m_density_field_id,
            "flattened density field");
    m_flux_field_id = field::addToStorage<FluxField>(
        m_lattice->get_blocks(), "flux field", FloatType{0}, field::fzyx,
        m_lattice->get_ghost_layers());
    m_flux_field_flattened_id =
        field::addFlattenedShallowCopyToStorage<FluxField>(
            m_lattice->get_blocks(), m_flux_field_id, "flattened flux field");

    m_continuity = std::make_unique<ContinuityKernel>(
        m_flux_field_flattened_id, m_density_field_flattened_id);

    // Init boundary related stuff
    m_flag_field_density_id = field::addFlagFieldToStorage<FlagField>(
        m_lattice->get_blocks(), "flag field density",
        m_lattice->get_ghost_layers());
    reset_density_boundary_handling();

    m_flag_field_flux_id = field::addFlagFieldToStorage<FlagField>(
        m_lattice->get_blocks(), "flag field flux",
        m_lattice->get_ghost_layers());
    reset_flux_boundary_handling();

    m_full_communication =
        std::make_shared<FullCommunicator>(m_lattice->get_blocks());
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<DensityField>>(
            m_density_field_id));

    // Synchronize ghost layers
    (*m_full_communication)();
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
  [[nodiscard]] bool is_double_precision() const noexcept override {
    return std::is_same<double, FloatType>::value;
  }

  void set_diffusion(double diffusion) noexcept override {
    m_diffusion = FloatType_c(diffusion);
  }
  void set_kT(double kT) noexcept override { m_kT = FloatType_c(kT); }
  void set_valency(double valency) noexcept override {
    m_valency = FloatType_c(valency);
  }
  void set_advection(bool advection) noexcept override {
    m_advection = advection;
  }
  void set_friction_coupling(bool friction_coupling) noexcept override {
    m_friction_coupling = friction_coupling;
  }
  void set_ext_efield(const Utils::Vector3d &field) noexcept override {
    m_ext_efield = field;
  }

  void ghost_communication() override { (*m_full_communication)(); }

private:
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
    auto kernel = DiffusiveFluxKernel(m_flux_field_flattened_id,
                                      m_density_field_flattened_id,
                                      FloatType_c(get_diffusion()));

    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  void kernel_advection(const std::size_t &velocity_id) {
    auto kernel =
        AdvectiveFluxKernel(m_flux_field_flattened_id, m_density_field_id,
                            BlockDataID(velocity_id));
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  void kernel_friction_coupling(const std::size_t &force_id) {
    auto kernel = FrictionCouplingKernel(
        BlockDataID(force_id), m_flux_field_flattened_id,
        FloatType_c(get_diffusion()), FloatType_c(get_kT()));
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  void kernel_diffusion_electrostatic(const std::size_t &potential_id) {
    auto const ext_field = get_ext_efield();
    auto kernel = DiffusiveFluxKernelElectrostatic(
        m_flux_field_flattened_id, BlockDataID(potential_id),
        m_density_field_flattened_id, FloatType_c(get_diffusion()),
        FloatType_c(ext_field[0]), FloatType_c(ext_field[1]),
        FloatType_c(ext_field[2]), FloatType_c(get_kT()),
        FloatType_c(get_valency()));
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
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
  void integrate(const std::size_t &potential_id,
                 const std::size_t &velocity_id,
                 const std::size_t &force_id) override {

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
      kernel_friction_coupling(force_id);
    }

    if (get_advection()) {
      if (velocity_id == walberla::BlockDataID{}) {
        throw std::runtime_error("Walberla EK: advection enabled but no "
                                 "velocity field accessible. velocity_id is " +
                                 std::to_string(velocity_id) +
                                 ". Hint: LB may be inactive.");
      }
      kernel_advection(velocity_id);
    }
    kernel_continuity();
    // is not necessary when reactions are done
    //    ghost_communication();

    // is this the expected behavior when reactions are included?
    kernel_boundary_density();

    // Handle VTK writers
    integrate_vtk_writers();
  }

  [[nodiscard]] std::size_t get_density_id() const noexcept override {
    static_assert(std::is_same_v<std::size_t, walberla::uint_t>);
    return static_cast<std::size_t>(m_density_field_id);
  }

  // Density
  bool set_node_density(const Utils::Vector3i &node, double density) override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc)
      return false;

    auto density_field =
        (*bc).block->template getData<DensityField>(m_density_field_id);
    density_field->get((*bc).cell) = FloatType_c(density);

    return true;
  }

  [[nodiscard]] boost::optional<double>
  get_node_density(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(get_lattice(), node, false);

    if (!bc)
      return {boost::none};

    auto density_field =
        (*bc).block->template getData<DensityField>(m_density_field_id);

    return {double_c(density_field->get((*bc).cell))};
  }

  void clear_flux_boundaries() override { reset_flux_boundary_handling(); }
  void clear_density_boundaries() override {
    reset_density_boundary_handling();
  }

  bool set_node_flux_boundary(const Utils::Vector3i &node,
                              const Utils::Vector3d &flux) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_flux->set_node_value_at_boundary(node, flux, *bc);

    return true;
  }

  [[nodiscard]] boost::optional<Utils::Vector3d>
  get_node_flux_at_boundary(const Utils::Vector3i &node) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc or !m_boundary_flux->node_is_boundary(node))
      return {boost::none};

    return {m_boundary_flux->get_node_value_at_boundary(node)};
  }

  bool remove_node_from_flux_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_flux->remove_node_from_boundary(node, *bc);

    return true;
  }

  bool set_node_density_boundary(const Utils::Vector3i &node,
                                 double density) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_density->set_node_value_at_boundary(node, FloatType_c(density),
                                                   *bc);

    return true;
  }

  [[nodiscard]] boost::optional<double>
  get_node_density_at_boundary(const Utils::Vector3i &node) const override {
    auto const bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc or !m_boundary_density->node_is_boundary(node))
      return {boost::none};

    return {double_c(m_boundary_density->get_node_value_at_boundary(node))};
  }

  bool remove_node_from_density_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_density->remove_node_from_boundary(node, *bc);

    return true;
  }

  [[nodiscard]] boost::optional<bool>
  get_node_is_flux_boundary(const Utils::Vector3i &node,
                            bool consider_ghosts) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    return {m_boundary_flux->node_is_boundary(node)};
  }

  [[nodiscard]] boost::optional<bool>
  get_node_is_density_boundary(const Utils::Vector3i &node,
                               bool consider_ghosts) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    return {m_boundary_density->node_is_boundary(node)};
  }

  [[nodiscard]] boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    return {m_boundary_density->node_is_boundary(node) or
            m_boundary_flux->node_is_boundary(node)};
  }

  void update_flux_boundary_from_shape(
      const std::vector<int> &raster_flat,
      const std::vector<double> &data_flat) override {
    auto const grid_size = get_lattice().get_grid_dimensions();
    auto const data = fill_3D_vector_array(data_flat, grid_size);
    auto const field_getter = [field_id = m_flux_field_id](auto &block) {
      return block->template getData<FluxField>(field_id);
    };
    set_boundary_from_grid(*m_boundary_flux, field_getter, get_lattice(),
                           raster_flat, data);
    reallocate_flux_boundary_field();
  }

  void update_density_boundary_from_shape(
      const std::vector<int> &raster_flat,
      const std::vector<double> &data_flat) override {
    auto const grid_size = get_lattice().get_grid_dimensions();
    auto const data = fill_3D_scalar_array(data_flat, grid_size);
    auto const field_getter = [field_id = m_density_field_id](auto &block) {
      return block->template getData<DensityField>(field_id);
    };
    set_boundary_from_grid(*m_boundary_density, field_getter, get_lattice(),
                           raster_flat, data);
    reallocate_density_boundary_field();
  }

  void reallocate_flux_boundary_field() { m_boundary_flux->boundary_update(); }

  void reallocate_density_boundary_field() {
    m_boundary_density->boundary_update();
  }

  [[nodiscard]] LatticeWalberla const &get_lattice() const noexcept override {
    return *m_lattice;
  }

  void register_vtk_field_filters(walberla::vtk::VTKOutput &vtk_obj) override {
    field::FlagFieldCellFilter<FlagField> fluid_filter(m_flag_field_density_id);
    fluid_filter.addFlag(Boundary_flag);
    vtk_obj.addCellExclusionFilter(fluid_filter);
  }

  void register_vtk_field_writers(walberla::vtk::VTKOutput &vtk_obj,
                                  int flag_observables) override {
    if (flag_observables & static_cast<int>(EKOutputVTK::density)) {
      vtk_obj.addCellDataWriter(
          make_shared<field::VTKWriter<DensityField, float>>(m_density_field_id,
                                                             "density"));
    }
  }

  ~EKinWalberlaImpl() override = default;
};

} // namespace walberla
