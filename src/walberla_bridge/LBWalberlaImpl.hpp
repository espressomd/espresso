/*
 * Copyright (C) 2020 The ESPResSo project
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
#ifndef LB_WALBERLA_H
#define LB_WALBERLA_H

/**
 * @file
 * @ref walberla::LBWalberlaImpl implements the interface of the LB
 * waLBerla bridge using sweeps generated by lbmpy
 * (see <tt>maintainer/walberla_kernels</tt>).
 */

#include <blockforest/Initialization.h>
#include <blockforest/StructuredBlockForest.h>
#include <blockforest/communication/UniformBufferedScheme.h>
#include <domain_decomposition/IBlock.h>
#include <field/GhostLayerField.h>
#include <field/vtk/FlagFieldCellFilter.h>
#include <field/vtk/VTKWriter.h>

#include <domain_decomposition/SharedSweep.h>

#include <field/AddToStorage.h>
#include <field/FlagField.h>
#include <field/communication/PackInfo.h>
#include <lbm/communication/PdfFieldPackInfo.h>
#include <lbm/field/AddToStorage.h>
#include <lbm/field/PdfField.h>
#include <lbm/sweeps/CellwiseSweep.h>

#include <stencil/D3Q19.h>
#include <stencil/D3Q27.h>

#include "BlockAndCell.hpp"
#include "BoundaryHandling.hpp"
#include "LBWalberlaBase.hpp"
#include "LatticeWalberla.hpp"
#include "LeesEdwardsPack.hpp"
#include "LeesEdwardsUpdate.hpp"
#include "ResetForce.hpp"
#include "generated_kernels/InitialPDFsSetterDoublePrecision.h"
#include "generated_kernels/InitialPDFsSetterSinglePrecision.h"
#include "generated_kernels/StreamSweepDoublePrecision.h"
#include "generated_kernels/StreamSweepSinglePrecision.h"
#include "generated_kernels/macroscopic_values_accessors_double_precision.h"
#include "generated_kernels/macroscopic_values_accessors_single_precision.h"
#include "vtk_writers.hpp"
#include "walberla_utils.hpp"

#ifdef __AVX2__
#include "generated_kernels/CollideSweepDoublePrecisionAVX.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalizedAVX.h"
#include "generated_kernels/CollideSweepLeesEdwardsAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalizedAVX.h"
#else
#include "generated_kernels/CollideSweepDoublePrecision.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalized.h"
#include "generated_kernels/CollideSweepLeesEdwards.h"
#include "generated_kernels/CollideSweepSinglePrecision.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalized.h"
#endif

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>
#include <utils/math/make_lin_space.hpp>

#include <boost/multi_array.hpp>
#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/variant.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace walberla {

namespace detail {
template <typename FT = double> struct KernelTrait {
#ifdef __AVX2__
  using ThermalizedCollisionModel =
      pystencils::CollideSweepDoublePrecisionThermalizedAVX;
  using UnthermalizedCollisionModel =
      pystencils::CollideSweepDoublePrecisionAVX;
  using LeesEdwardsCollisionModel = pystencils::CollideSweepLeesEdwardsAVX;
#else
  using ThermalizedCollisionModel =
      pystencils::CollideSweepDoublePrecisionThermalized;
  using UnthermalizedCollisionModel = pystencils::CollideSweepDoublePrecision;
  using LeesEdwardsCollisionModel = pystencils::CollideSweepLeesEdwards;
#endif
  using StreamSweep = pystencils::StreamSweepDoublePrecision;
  using InitialPDFsSetter = pystencils::InitialPDFsSetterDoublePrecision;
};
template <> struct KernelTrait<float> {
#ifdef __AVX2__
  using ThermalizedCollisionModel =
      pystencils::CollideSweepSinglePrecisionThermalizedAVX;
  using UnthermalizedCollisionModel =
      pystencils::CollideSweepSinglePrecisionAVX;
  using LeesEdwardsCollisionModel = pystencils::CollideSweepLeesEdwardsAVX;
#else
  using ThermalizedCollisionModel =
      pystencils::CollideSweepSinglePrecisionThermalized;
  using UnthermalizedCollisionModel = pystencils::CollideSweepSinglePrecision;
  using LeesEdwardsCollisionModel = pystencils::CollideSweepLeesEdwards;
#endif
  using StreamSweep = pystencils::StreamSweepSinglePrecision;
  using InitialPDFsSetter = pystencils::InitialPDFsSetterSinglePrecision;
};
} // namespace detail

/** Class that runs and controls the LB on WaLBerla
 */
template <typename FloatType = double>
class LBWalberlaImpl : public LBWalberlaBase {
  template <typename T> inline FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }
  using LeesEdwardsCollisionModel =
      typename detail::KernelTrait<FloatType>::LeesEdwardsCollisionModel;
  using UnthermalizedCollisionModel =
      typename detail::KernelTrait<FloatType>::UnthermalizedCollisionModel;
  using ThermalizedCollisionModel =
      typename detail::KernelTrait<FloatType>::ThermalizedCollisionModel;
  using StreamSweep = typename detail::KernelTrait<FloatType>::StreamSweep;
  using InitialPDFsSetter =
      typename detail::KernelTrait<FloatType>::InitialPDFsSetter;
  using BoundaryModel = BoundaryHandling<FloatType>;

protected:
  using CollisionModel =
      boost::variant<UnthermalizedCollisionModel, ThermalizedCollisionModel,
                     LeesEdwardsCollisionModel>;

private:
  class : public boost::static_visitor<> {
  public:
    void operator()(UnthermalizedCollisionModel &cm, IBlock *b) { cm(b); }

    void operator()(ThermalizedCollisionModel &cm, IBlock *b) { cm(b); }

    void operator()(LeesEdwardsCollisionModel &cm, IBlock *b) {
      //cm.shear_velocity_ = m_lees_edwards_sweep->get_shear_velocity();
      cm.points_up_ = m_lees_edwards_sweep->points_up(b);
      cm.points_down_ = m_lees_edwards_sweep->points_down(b);
      cm(b);
    }

    std::shared_ptr<LeesEdwardsUpdate> m_lees_edwards_sweep;

  } run_collide_sweep;

  FloatType shear_mode_relaxation_rate() const {
    return FloatType{2} / (FloatType{6} * m_viscosity + FloatType{1});
  }

  FloatType odd_mode_relaxation_rate(
      FloatType shear_relaxation,
      FloatType magic_number = FloatType{3} / FloatType{16}) const {
    return (FloatType{4} - FloatType{2} * shear_relaxation) /
           (FloatType{4} * magic_number * shear_relaxation + FloatType{2} -
            shear_relaxation);
  }

  void reset_boundary_handling() {
    auto const &blocks = lattice().get_blocks();
    m_boundary = std::make_shared<BoundaryModel>(blocks, m_pdf_field_id,
                                                 m_flag_field_id);
  }

public:
  // Type definitions
  typedef stencil::D3Q19 Stencil;
  using VectorField = GhostLayerField<FloatType, 3u>;
  using FlagField = typename BoundaryModel::FlagField;
  using PdfField = GhostLayerField<FloatType, Stencil::Size>;

private:
  FloatType getDensity(const BlockAndCell &bc) const {
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    return lbm::accessor::Density::get(*pdf_field, bc.cell.x(), bc.cell.y(),
                                       bc.cell.z());
  }

  FloatType getDensityAndVelocity(const BlockAndCell &bc,
                                  Vector3<FloatType> &velocity) const {
    auto const pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    auto const force_field =
        bc.block->template getData<VectorField>(m_last_applied_force_field_id);
    return getDensityAndVelocity(pdf_field, force_field, bc.cell.x(),
                                 bc.cell.y(), bc.cell.z(), velocity);
  }

  FloatType getDensityAndVelocity(const PdfField *pdf_field,
                                  const VectorField *force_field,
                                  const cell_idx_t x, const cell_idx_t y,
                                  const cell_idx_t z,
                                  Vector3<FloatType> &velocity) const {
    auto const rho = lbm::accessor::DensityAndMomentumDensity::get(
        velocity, *force_field, *pdf_field, x, y, z);
    auto const invRho = FloatType{1} / rho;
    velocity *= invRho;
    return rho;
  }

  void setDensityAndVelocity(const BlockAndCell &bc,
                             Vector3<FloatType> const &velocity,
                             FloatType rho) {
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    auto force_field =
        bc.block->template getData<VectorField>(m_last_applied_force_field_id);
    lbm::accessor::DensityAndVelocity::set(*pdf_field, bc.cell.x(), bc.cell.y(),
                                           bc.cell.z(), *force_field, velocity,
                                           rho);
  }

  Matrix3<FloatType> getPressureTensor(const BlockAndCell &bc) const {
    Matrix3<FloatType> pressureTensor;
    auto pdf_field = bc.block->template getData<PdfField>(m_pdf_field_id);
    lbm::accessor::PressureTensor::get(pressureTensor, *pdf_field, bc.cell.x(),
                                       bc.cell.y(), bc.cell.z());
    return pressureTensor;
  }

  class interpolation_illegal_access : public std::runtime_error {
  public:
    explicit interpolation_illegal_access(std::string const &field,
                                          Utils::Vector3d const &pos,
                                          std::array<int, 3> const &node,
                                          double weight)
        : std::runtime_error("Access to LB " + field + " field failed") {
      std::cerr << "pos [" << pos << "], "
                << "node [" << Utils::Vector3i(node) << "], "
                << "weight " << weight << "\n";
    }
  };

  class vtk_runtime_error : public std::runtime_error {
  public:
    explicit vtk_runtime_error(std::string const &vtk_uid,
                               std::string const &reason)
        : std::runtime_error("VTKOutput object '" + vtk_uid + "' " + reason) {}
  };

protected:
  /** VTK writers that are executed automatically */
  std::map<std::string, std::shared_ptr<VTKHandle>> m_vtk_auto;
  /** VTK writers that are executed manually */
  std::map<std::string, std::shared_ptr<VTKHandle>> m_vtk_manual;

  // Member variables
  FloatType m_viscosity;
  FloatType m_density;
  FloatType m_kT;

  // Block data access handles
  BlockDataID m_pdf_field_id;
  BlockDataID m_pdf_tmp_field_id;
  BlockDataID m_flag_field_id;

  BlockDataID m_last_applied_force_field_id;
  BlockDataID m_force_to_be_applied_id;

  BlockDataID m_velocity_field_id;

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;
  using PDFStreamingCommunicator =
      blockforest::communication::UniformBufferedScheme<
          typename stencil::D3Q19>;
  std::shared_ptr<PDFStreamingCommunicator> m_pdf_streaming_communication;

  // ResetForce sweep + external force handling
  std::shared_ptr<ResetForce<PdfField, VectorField>> m_reset_force;

  // Stream sweep
  std::shared_ptr<StreamSweep> m_stream;

  // Collision sweep
  std::shared_ptr<CollisionModel> m_collision_model;

  // boundaries
  std::shared_ptr<BoundaryModel> m_boundary;

  // lattice
  std::shared_ptr<LatticeWalberla> m_lattice;

  // Lees-Edwards sweep
  std::shared_ptr<LeesEdwardsUpdate> m_lees_edwards_sweep;

  std::size_t stencil_size() const override {
    return static_cast<std::size_t>(Stencil::Size);
  }

public:
  LBWalberlaImpl(std::shared_ptr<LatticeWalberla> const &lattice,
                 double viscosity, double density)
      : m_viscosity(FloatType_c(viscosity)), m_density(FloatType_c(density)),
        m_kT(FloatType{0}), m_lattice(lattice) {

    auto const &blocks = m_lattice->get_blocks();
    auto const n_ghost_layers = m_lattice->get_ghost_layers();
    if (n_ghost_layers == 0)
      throw std::runtime_error("At least one ghost layer must be used");

    // Init and register fields
    m_pdf_field_id = field::addToStorage<PdfField>(blocks, "pdfs", FloatType{0},
                                                   field::fzyx, n_ghost_layers);
    m_pdf_tmp_field_id = field::addToStorage<PdfField>(
        blocks, "pdfs_tmp", FloatType{0}, field::fzyx, n_ghost_layers);
    m_last_applied_force_field_id = field::addToStorage<VectorField>(
        blocks, "force field", FloatType{0}, field::fzyx, n_ghost_layers);
    m_force_to_be_applied_id = field::addToStorage<VectorField>(
        blocks, "force field", FloatType{0}, field::fzyx, n_ghost_layers);
    m_velocity_field_id = field::addToStorage<VectorField>(
        blocks, "velocity field", FloatType{0}, field::fzyx, n_ghost_layers);

    // Init and register pdf field
    auto pdf_setter =
        InitialPDFsSetter(m_force_to_be_applied_id, m_pdf_field_id,
                          m_velocity_field_id, m_density);
    for (auto b = blocks->begin(); b != blocks->end(); ++b) {
      pdf_setter(&*b);
    }

    // Init and register flag field (fluid/boundary)
    m_flag_field_id = field::addFlagFieldToStorage<FlagField>(
        blocks, "flag field", n_ghost_layers);
    // Init boundary sweep
    reset_boundary_handling();

    // Set up the communication and register fields
    m_pdf_streaming_communication =
        std::make_shared<PDFStreamingCommunicator>(blocks);
    m_pdf_streaming_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PdfField>>(
            m_pdf_field_id, n_ghost_layers));
    m_pdf_streaming_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_last_applied_force_field_id, n_ghost_layers));
    m_pdf_streaming_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<FlagField>>(
            m_flag_field_id, n_ghost_layers));

    m_full_communication = std::make_shared<FullCommunicator>(blocks);
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PdfField>>(
            m_pdf_field_id, n_ghost_layers));
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_last_applied_force_field_id, n_ghost_layers));
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<VectorField>>(
            m_velocity_field_id, n_ghost_layers));
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<FlagField>>(
            m_flag_field_id, n_ghost_layers));

    // Instance the sweep responsible for force double buffering and
    // external forces
    m_reset_force = std::make_shared<ResetForce<PdfField, VectorField>>(
        m_last_applied_force_field_id, m_force_to_be_applied_id);

    // Prepare LB sweeps
    // Note: For now, combined collide-stream sweeps cannot be used,
    // because the collide-push variant is not supported by lbmpy.
    // The following functors are individual in-place collide and stream steps
    m_stream = std::make_shared<StreamSweep>(
        m_last_applied_force_field_id, m_pdf_field_id, m_velocity_field_id);
    set_collision_model();

    // Synchronize ghost layers
    (*m_full_communication)();
  }

  void integrate() override {
    auto const &blocks = lattice().get_blocks();
    // Reset force fields
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_reset_force)(&*b);
    // Handle boundaries
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_boundary)(&*b);
    // LB stream
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      (*m_stream)(&*b);
    // Refresh ghost layers
    (*m_full_communication)();
    // LB collide
    for (auto b = blocks->begin(); b != blocks->end(); ++b)
      boost::apply_visitor(run_collide_sweep, *m_collision_model,
                           boost::variant<IBlock *>(&*b));
    if (auto *cm = boost::get<ThermalizedCollisionModel>(&*m_collision_model)) {
      cm->time_step_++;
    }

    // Handle VTK writers
    for (auto it = m_vtk_auto.begin(); it != m_vtk_auto.end(); ++it) {
      auto &vtk_handle = it->second;
      if (vtk_handle->enabled) {
        vtk::writeFiles(vtk_handle->ptr)();
        vtk_handle->execution_count++;
      }
    }
  }

  void ghost_communication() override { (*m_full_communication)(); }

  void set_collision_model() override {
    auto const omega = shear_mode_relaxation_rate();
    auto const omega_odd = odd_mode_relaxation_rate(omega);
    auto obj = UnthermalizedCollisionModel(m_last_applied_force_field_id,
                                           m_pdf_field_id, omega, omega,
                                           omega_odd, omega);
    m_collision_model = std::make_shared<CollisionModel>(std::move(obj));
  }

  void set_collision_model(double kT, unsigned int seed) override {
    auto const omega = shear_mode_relaxation_rate();
    auto const omega_odd = odd_mode_relaxation_rate(omega);
    m_kT = FloatType_c(kT);
    auto obj = ThermalizedCollisionModel(
        m_last_applied_force_field_id, m_pdf_field_id, uint32_t(0u),
        uint32_t(0u), uint32_t(0u), m_kT, omega, omega, omega_odd, omega, seed,
        0u);
    obj.block_offset_generator =
        [this](IBlock *const block, uint32_t &block_offset_0,
               uint32_t &block_offset_1, uint32_t &block_offset_2) {
          auto const &blocks = lattice().get_blocks();
          block_offset_0 = blocks->getBlockCellBB(*block).xMin();
          block_offset_1 = blocks->getBlockCellBB(*block).yMin();
          block_offset_2 = blocks->getBlockCellBB(*block).zMin();
        };
    m_collision_model = std::make_shared<CollisionModel>(std::move(obj));
  }

  void set_collision_model(LeesEdwardsPack &&lees_edwards_pack) override {
    if (m_kT != 0.) {
      throw std::runtime_error(
          "Lees-Edwards LB doesn't support thermalization");
    }
    auto const shear_plane_size =
        Utils::product(lattice().get_grid_dimensions()) /
        lattice().get_grid_dimensions()[lees_edwards_pack.shear_plane_normal];
    m_lees_edwards_sweep = std::make_shared<LeesEdwardsUpdate>(
        lattice().get_blocks(), m_pdf_field_id, m_pdf_tmp_field_id,
        lattice().get_ghost_layers(), std::move(lees_edwards_pack));
    run_collide_sweep.m_lees_edwards_sweep = m_lees_edwards_sweep;
    auto const omega = shear_mode_relaxation_rate();
    auto const omega_odd = odd_mode_relaxation_rate(omega);
    // a few values are initialized to 0 or false, will be updated later
    auto obj = LeesEdwardsCollisionModel(
        m_last_applied_force_field_id, m_pdf_field_id, m_velocity_field_id,
        omega, omega, omega_odd, omega, false, false);
    m_collision_model = std::make_shared<CollisionModel>(std::move(obj));
    //auto *cm = boost::get<LeesEdwardsCollisionModel>(&*m_collision_model);
    //cm->grid_size_ = int64_t(shear_plane_size);
  }

  void set_viscosity(double viscosity) override {
    m_viscosity = FloatType_c(viscosity);
  }

  double get_viscosity() const override { return m_viscosity; }

  double get_density() const override { return m_density; }

  // Velocity
  boost::optional<Utils::Vector3d>
  get_node_velocity(const Utils::Vector3i &node,
                    bool consider_ghosts = false) const override {
    auto const is_boundary = get_node_is_boundary(node, consider_ghosts);
    if (is_boundary)    // is info available locally
      if (*is_boundary) // is the node a boundary
        return get_node_velocity_at_boundary(node);
    auto const bc = get_block_and_cell(lattice(), node, consider_ghosts);
    if (!bc)
      return {};
    auto const &vel_field =
        (*bc).block->template getData<VectorField>(m_velocity_field_id);
    return Utils::Vector3d{double_c(vel_field->get((*bc).cell, uint_t(0u))),
                           double_c(vel_field->get((*bc).cell, uint_t(1u))),
                           double_c(vel_field->get((*bc).cell, uint_t(2u)))};
  }
  bool set_node_velocity(const Utils::Vector3i &node,
                         const Utils::Vector3d &v) override {
    auto bc = get_block_and_cell(lattice(), node, false);
    if (!bc)
      return false;
    // We have to set both, the pdf and the stored velocity field
    auto const density = getDensity(*bc);
    auto const vel = to_vector3<FloatType>(v);
    setDensityAndVelocity(*bc, vel, density);
    auto vel_field =
        (*bc).block->template getData<VectorField>(m_velocity_field_id);
    for (uint_t f = 0u; f < 3u; ++f) {
      vel_field->get((*bc).cell, f) = FloatType_c(v[f]);
    }

    return true;
  }
  boost::optional<Utils::Vector3d>
  get_velocity_at_pos(const Utils::Vector3d &pos,
                      bool consider_points_in_halo = false) const override {
    if (!consider_points_in_halo and !m_lattice->pos_in_local_domain(pos))
      return {};
    if (consider_points_in_halo and !m_lattice->pos_in_local_halo(pos))
      return {};
    Utils::Vector3d v{0.0, 0.0, 0.0};
    interpolate_bspline_at_pos(
        pos, [this, &v, pos](const std::array<int, 3> node, double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0) {
            auto res = get_node_velocity(
                Utils::Vector3i{{node[0], node[1], node[2]}}, true);
            if (!res) {
              throw interpolation_illegal_access("velocity", pos, node, weight);
            }
            v += *res * weight;
          }
        });
    return {v};
  }

  boost::optional<double> get_interpolated_density_at_pos(
      const Utils::Vector3d &pos,
      bool consider_points_in_halo = false) const override {
    if (!consider_points_in_halo and !m_lattice->pos_in_local_domain(pos))
      return {};
    if (consider_points_in_halo and !m_lattice->pos_in_local_halo(pos))
      return {};
    double dens = 0.0;
    interpolate_bspline_at_pos(
        pos, [this, &dens, pos](const std::array<int, 3> node, double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0) {
            auto const res = get_node_density(Utils::Vector3i(node));
            if (!res) {
              throw interpolation_illegal_access("density", pos, node, weight);
            }
            dens += *res * weight;
          }
        });
    return {dens};
  }

  // Local force
  bool add_force_at_pos(const Utils::Vector3d &pos,
                        const Utils::Vector3d &force) override {
    if (!m_lattice->pos_in_local_halo(pos))
      return false;
    auto force_at_node = [this, force](const std::array<int, 3> node,
                                       double weight) {
      auto const bc = get_block_and_cell(lattice(), to_vector3i(node), true);
      if (bc) {
        auto force_field = (*bc).block->template getData<VectorField>(
            m_force_to_be_applied_id);
        for (int i : {0, 1, 2})
          force_field->get((*bc).cell, i) += FloatType_c(force[i] * weight);
      }
    };
    interpolate_bspline_at_pos(pos, force_at_node);
    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_force_to_be_applied(const Utils::Vector3i &node) const override {
    auto const bc = get_block_and_cell(lattice(), node, true);
    if (!bc)
      return {};

    auto const &force_field =
        (*bc).block->template getData<VectorField>(m_force_to_be_applied_id);
    return Utils::Vector3d{double_c(force_field->get((*bc).cell, uint_t(0u))),
                           double_c(force_field->get((*bc).cell, uint_t(1u))),
                           double_c(force_field->get((*bc).cell, uint_t(2u)))};
  }

  bool set_node_last_applied_force(Utils::Vector3i const &node,
                                   Utils::Vector3d const &force) override {
    auto bc = get_block_and_cell(lattice(), node, false);
    if (!bc)
      return false;

    auto force_field = (*bc).block->template getData<VectorField>(
        m_last_applied_force_field_id);
    for (uint_t f = 0u; f < 3u; ++f) {
      force_field->get((*bc).cell, f) = FloatType_c(force[f]);
    }

    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_last_applied_force(const Utils::Vector3i &node,
                              bool consider_ghosts = false) const override {
    auto const bc = get_block_and_cell(lattice(), node, consider_ghosts);
    if (!bc)
      return {};

    auto const force_field = (*bc).block->template getData<VectorField>(
        m_last_applied_force_field_id);
    return Utils::Vector3d{double_c(force_field->get((*bc).cell, uint_t(0u))),
                           double_c(force_field->get((*bc).cell, uint_t(1u))),
                           double_c(force_field->get((*bc).cell, uint_t(2u)))};
  }

  // Population
  bool set_node_pop(Utils::Vector3i const &node,
                    std::vector<double> const &population) override {
    auto bc = get_block_and_cell(lattice(), node, false);
    if (!bc)
      return false;

    auto pdf_field = (*bc).block->template getData<PdfField>(m_pdf_field_id);
    assert(population.size() == Stencil::Size);
    for (uint_t f = 0u; f < Stencil::Size; ++f) {
      pdf_field->get((*bc).cell, f) = FloatType_c(population[f]);
    }

    return true;
  }

  boost::optional<std::vector<double>>
  get_node_pop(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(lattice(), node, false);
    if (!bc)
      return {boost::none};

    auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
    std::vector<double> population(Stencil::Size);
    for (uint_t f = 0u; f < Stencil::Size; ++f) {
      population[f] = double_c(pdf_field->get((*bc).cell, f));
    }

    return {population};
  }

  // Density
  bool set_node_density(const Utils::Vector3i &node, double density) override {
    auto bc = get_block_and_cell(lattice(), node, false);
    if (!bc)
      return false;

    Vector3<FloatType> vel;
    getDensityAndVelocity(*bc, vel);
    setDensityAndVelocity(*bc, vel, FloatType_c(density));

    return true;
  }

  boost::optional<double>
  get_node_density(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(lattice(), node, false);
    if (!bc)
      return {boost::none};

    auto const density = getDensity(*bc);
    return {double_c(density)};
  }

  boost::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(lattice(), node, true);
    if (!bc or !m_boundary->node_is_boundary(*bc))
      return {boost::none};

    return {m_boundary->get_node_velocity_at_boundary(node)};
  }

  bool set_node_velocity_at_boundary(const Utils::Vector3i &node,
                                     const Utils::Vector3d &v,
                                     bool reallocate) override {
    auto bc = get_block_and_cell(lattice(), node, true);
    if (!bc)
      return false;

    m_boundary->set_node_velocity_at_boundary(node, v, *bc);
    if (reallocate) {
      m_boundary->ubb_update();
    }

    return true;
  }

  boost::optional<Utils::Vector3d>
  get_node_boundary_force(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(lattice(), node, true); // including ghosts
    if (!bc or !m_boundary->node_is_boundary(*bc))
      return {boost::none};

    return get_node_last_applied_force(node, true);
  }

  bool remove_node_from_boundary(const Utils::Vector3i &node,
                                 bool reallocate) override {
    auto bc = get_block_and_cell(lattice(), node, true);
    if (!bc)
      return false;

    m_boundary->remove_node_from_boundary(node, *bc);
    if (reallocate) {
      m_boundary->ubb_update();
    }

    return true;
  }

  boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    return {m_boundary->node_is_boundary(*bc)};
  }

  void reallocate_ubb_field() override { m_boundary->ubb_update(); }

  void clear_boundaries() override { reset_boundary_handling(); }

  /** @brief Update boundary conditions from a rasterized shape. */
  void update_boundary_from_shape(
      std::vector<int> const &raster_flat,
      std::vector<double> const &slip_velocity_flat) override {
    // reshape grids
    auto const grid_size = lattice().get_grid_dimensions();
    auto const n_grid_points = Utils::product(grid_size);
    assert(raster_flat.size() == n_grid_points);
    assert(slip_velocity_flat.size() == 3 * n_grid_points or
           slip_velocity_flat.size() == 3);
    std::vector<Utils::Vector3d> slip_velocity_vectors;
    {
      auto const vel_begin = std::begin(slip_velocity_flat);
      auto const vel_end = std::end(slip_velocity_flat);
      if (slip_velocity_flat.size() == 3) {
        auto const uniform_slip_velocity = Utils::Vector3d(vel_begin, vel_end);
        slip_velocity_vectors.assign(n_grid_points, uniform_slip_velocity);
      } else {
        slip_velocity_vectors.reserve(n_grid_points);
        for (auto it = vel_begin; it < vel_end; it += 3) {
          slip_velocity_vectors.emplace_back(Utils::Vector3d(it, it + 3));
        }
      }
    }
    boost::const_multi_array_ref<Utils::Vector3d, 3> slip_velocity(
        slip_velocity_vectors.data(), grid_size);
    boost::const_multi_array_ref<int, 3> raster(raster_flat.data(), grid_size);

    auto const &blocks = lattice().get_blocks();
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
      // lattice constant is 1
      auto const left = block->getAABB().min();
      auto const off_i = int_c(left[0]);
      auto const off_j = int_c(left[1]);
      auto const off_k = int_c(left[2]);

      // Get field data which knows about the indices
      // In the loop, x,y,z are in block-local coordinates
      auto const n_ghost_layers = lattice().get_ghost_layers();
      auto const ghosts = static_cast<int>(n_ghost_layers);
      auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
      for (int i = -ghosts; i < int_c(pdf_field->xSize() + ghosts); ++i) {
        for (int j = -ghosts; j < int_c(pdf_field->ySize() + ghosts); ++j) {
          for (int k = -ghosts; k < int_c(pdf_field->zSize() + ghosts); ++k) {
            Utils::Vector3i const node{{off_i + i, off_j + j, off_k + k}};
            auto const idx = (node + grid_size) % grid_size;
            if (raster(idx)) {
              auto const bc = get_block_and_cell(lattice(), node, true);
              auto const &vel = slip_velocity(idx);
              m_boundary->set_node_velocity_at_boundary(node, vel, *bc);
            }
          }
        }
      }
    }
    reallocate_ubb_field();
  }

  /** @brief Update boundary conditions from a list of nodes. */
  void update_boundary_from_list(std::vector<int> const &nodes_flat,
                                 std::vector<double> const &vel_flat) override {
    // reshape grids
    auto const grid_size = lattice().get_grid_dimensions();
    auto const n_grid_points = Utils::product(grid_size);
    auto const n_ghost_layers = lattice().get_ghost_layers();
    assert(nodes_flat.size() == vel_flat.size());
    assert(nodes_flat.size() % 3u == 0);
    for (std::size_t i = 0; i < nodes_flat.size(); i += 3) {
      auto const node = Utils::Vector3i(&nodes_flat[i], &nodes_flat[i + 3]);
      auto const vel = Utils::Vector3d(&vel_flat[i], &vel_flat[i + 3]);
      auto const bc = get_block_and_cell(lattice(), node, true);
      if (bc) {
        m_boundary->set_node_velocity_at_boundary(node, vel, *bc);
      }
    }
    reallocate_ubb_field();
  }

  // Pressure tensor
  boost::optional<Utils::Vector6d>
  get_node_pressure_tensor(const Utils::Vector3i &node) const override {
    auto const n_ghost_layers = lattice().get_ghost_layers();
    auto bc = get_block_and_cell(lattice(), node, false);
    if (!bc)
      return {boost::none};
    return to_vector6d(getPressureTensor(*bc));
  }

  // Global momentum
  Utils::Vector3d get_momentum() const override {
    auto const &blocks = lattice().get_blocks();
    Vector3<FloatType> mom;
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
      auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
      auto force_field =
          block->template getData<VectorField>(m_last_applied_force_field_id);
      Vector3<FloatType> local_v;
      WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
        FloatType local_dens =
            getDensityAndVelocity(pdf_field, force_field, x, y, z, local_v);
        mom += local_dens * local_v;
      });
    }
    return to_vector3d(mom);
  }
  // Global external force
  void set_external_force(const Utils::Vector3d &ext_force) override {
    m_reset_force->set_ext_force(ext_force);
  }
  Utils::Vector3d get_external_force() const override {
    return m_reset_force->get_ext_force();
  }

  double get_kT() const override { return m_kT; }

  uint64_t get_rng_state() const override {
    auto const *cm = boost::get<ThermalizedCollisionModel>(&*m_collision_model);
    if (!cm)
      throw std::runtime_error("The LB does not use a random number generator");
    return cm->time_step_;
  }
  void set_rng_state(uint64_t counter) override {
    auto *cm = boost::get<ThermalizedCollisionModel>(&*m_collision_model);
    if (!cm)
      throw std::runtime_error("The LB does not use a random number generator");
    cm->time_step_ = counter;
  }

  LatticeWalberla const &lattice() const override { return *m_lattice; }

  std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>>
  node_indices_positions(bool include_ghosts = false) const override {
    int ghost_offset = 0;
    if (include_ghosts)
      ghost_offset = lattice().get_ghost_layers();
    auto const &blocks = lattice().get_blocks();
    std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>> res;
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
      auto left = block->getAABB().min();
      // Lattice constant is 1, node centers are offset by .5
      Utils::Vector3d pos_offset =
          to_vector3d(left) + Utils::Vector3d::broadcast(.5);

      // Lattice constant is 1, so cast left corner position to ints
      Utils::Vector3i index_offset =
          Utils::Vector3i{int(left[0]), int(left[1]), int(left[2])};

      // Get field data which knows about the indices
      // In the loop, x,y,z are in block-local coordinates
      auto pdf_field = block->template getData<PdfField>(m_pdf_field_id);
      for (int x = -ghost_offset; x < int(pdf_field->xSize()) + ghost_offset;
           x++) {
        for (int y = -ghost_offset; y < int(pdf_field->ySize()) + ghost_offset;
             y++) {
          for (int z = -ghost_offset;
               z < int(pdf_field->zSize()) + ghost_offset; z++) {
            res.push_back({index_offset + Utils::Vector3i{x, y, z},
                           pos_offset + Utils::Vector3d{double(x), double(y),
                                                        double(z)}});
          }
        }
      }
    }
    return res;
  }

  std::shared_ptr<VTKHandle> create_vtk(int delta_N, int initial_count,
                                        int flag_observables,
                                        std::string const &identifier,
                                        std::string const &base_folder,
                                        std::string const &prefix) override {
    // VTKOuput object must be unique
    std::stringstream unique_identifier;
    unique_identifier << base_folder << "/" << identifier;
    std::string const vtk_uid = unique_identifier.str();
    if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end() or
        m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
      throw vtk_runtime_error(vtk_uid, "already exists");
    }

    // instantiate VTKOutput object
    auto const &blocks = lattice().get_blocks();
    auto const write_freq = (delta_N) ? static_cast<unsigned int>(delta_N) : 1u;
    auto pdf_field_vtk = vtk::createVTKOutput_BlockData(
        blocks, identifier, uint_c(write_freq), uint_c(0), false, base_folder,
        prefix, true, true, true, true, uint_c(initial_count));
    field::FlagFieldCellFilter<FlagField> fluid_filter(m_flag_field_id);
    fluid_filter.addFlag(Boundary_flag);
    pdf_field_vtk->addCellExclusionFilter(fluid_filter);

    // add writers
    if (flag_observables & static_cast<int>(OutputVTK::density)) {
      pdf_field_vtk->addCellDataWriter(
          make_shared<lbm::DensityVTKWriter<LBWalberlaImpl, float>>(
              m_pdf_field_id, "DensityFromPDF"));
    }
    if (flag_observables & static_cast<int>(OutputVTK::velocity_vector)) {
      pdf_field_vtk->addCellDataWriter(
          make_shared<field::VTKWriter<VectorField, float>>(
              m_velocity_field_id, "VelocityFromVelocityField"));
    }
    if (flag_observables & static_cast<int>(OutputVTK::pressure_tensor)) {
      pdf_field_vtk->addCellDataWriter(
          make_shared<lbm::PressureTensorVTKWriter<LBWalberlaImpl, float>>(
              m_pdf_field_id, "PressureTensorFromPDF"));
    }

    // register object
    auto vtk_handle =
        std::make_shared<VTKHandle>(pdf_field_vtk, initial_count, true);
    if (delta_N) {
      m_vtk_auto[vtk_uid] = vtk_handle;
    } else {
      m_vtk_manual[vtk_uid] = vtk_handle;
    }
    return vtk_handle;
  }

  /** Manually call a VTK callback */
  void write_vtk(std::string const &vtk_uid) override {
    if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end()) {
      throw vtk_runtime_error(vtk_uid, "is an automatic observable");
    }
    if (m_vtk_manual.find(vtk_uid) == m_vtk_manual.end()) {
      throw vtk_runtime_error(vtk_uid, "doesn't exist");
    }
    auto &vtk_handle = m_vtk_manual[vtk_uid];
    vtk::writeFiles(vtk_handle->ptr)();
    vtk_handle->execution_count++;
  }

  /** Activate or deactivate a VTK callback */
  void switch_vtk(std::string const &vtk_uid, bool status) override {
    if (m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
      throw vtk_runtime_error(vtk_uid, "is a manual observable");
    }
    if (m_vtk_auto.find(vtk_uid) == m_vtk_auto.end()) {
      throw vtk_runtime_error(vtk_uid, "doesn't exist");
    }
    m_vtk_auto[vtk_uid]->enabled = status;
  }

  ~LBWalberlaImpl() override = default;
};
} // namespace walberla

#endif // LB_WALBERLA_H
