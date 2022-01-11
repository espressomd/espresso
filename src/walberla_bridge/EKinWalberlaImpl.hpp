#ifndef ESPRESSO_EKINWALBERLAIMPL_HPP
#define ESPRESSO_EKINWALBERLAIMPL_HPP

#include "blockforest/communication/UniformBufferedScheme.h"
#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/lattice_model/D3Q27.h"
#include "timeloop/SweepTimeloop.h"

#include "BlockAndCell.hpp"
#include "EKinWalberlaBase.hpp"
#include "LatticeWalberla.hpp"
#include "walberla_utils.hpp"

#include "generated_kernels/AdvectiveFluxKernel.h"
#include "generated_kernels/ContinuityKernel.h"
#include "generated_kernels/DiffusiveFluxKernel.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostatic.h"
#include "generated_kernels/Dirichlet.h"
#include "generated_kernels/FixedFlux.h"
#include "generated_kernels/FrictionCouplingKernel.h"

#include <boost/multi_array/multi_array_ref.hpp>
#include <memory>

#include "BoundaryHandling.hpp"

namespace walberla {

std::vector<Utils::Vector3d>
fill_3D_vector_array(const Utils::Vector3i &grid_size,
                     const std::vector<double> &vec_flat) {
  auto const n_grid_points = Utils::product(grid_size);
  assert(vec_flat.size() == 3 * n_grid_points or vec_flat.size() == 3);
  std::vector<Utils::Vector3d> output_vector;

  auto const vec_begin = std::begin(vec_flat);
  auto const vec_end = std::end(vec_flat);
  if (vec_flat.size() == 3) {
    auto const uniform_vector = Utils::Vector3d(vec_begin, vec_end);
    output_vector.assign(n_grid_points, uniform_vector);
  } else {
    output_vector.reserve(n_grid_points);
    for (auto it = vec_begin; it < vec_end; it += 3) {
      output_vector.emplace_back(Utils::Vector3d(it, it + 3));
    }
  }

  return output_vector;
}

std::vector<double> fill_3D_scalar_array(const Utils::Vector3i &grid_size,
                                         const std::vector<double> &vec_flat) {
  auto const n_grid_points = Utils::product(grid_size);
  assert(vec_flat.size() == n_grid_points or vec_flat.size() == 1);
  std::vector<double> output_vector;

  auto const vel_begin = std::begin(output_vector);
  auto const vel_end = std::end(output_vector);
  if (vec_flat.size() == 1) {
    auto const uniform_value = vec_flat[0];
    output_vector.assign(n_grid_points, uniform_value);
  } else {
    output_vector.assign(vel_begin, vel_end);
  }

  return output_vector;
}

/** Class that runs and controls the LB on WaLBerla
 */
template <size_t FluxCount = 13, typename FloatType = double>
class EKinWalberlaImpl : public EKinWalberlaBase<FloatType> {
protected:
  // Type definitions
  using FluxField = GhostLayerField<FloatType, FluxCount>;
  using FlagField = walberla::FlagField<walberla::uint8_t>;
  using DensityField = GhostLayerField<FloatType, 1>;

  using EKinWalberlaBase<FloatType>::get_diffusion;
  using EKinWalberlaBase<FloatType>::get_kT;
  using EKinWalberlaBase<FloatType>::get_valency;
  using EKinWalberlaBase<FloatType>::get_ext_efield;
  using EKinWalberlaBase<FloatType>::get_advection;
  using EKinWalberlaBase<FloatType>::get_friction_coupling;

  // TODO: boundary handling for Density!
  using BoundaryModelDensity =
      BoundaryHandling<FloatType, pystencils::Dirichlet>;
  using BoundaryModelFlux =
      BoundaryHandling<Vector3<FloatType>, pystencils::FixedFlux>;

protected:
  /** VTK writers that are executed automatically */
  std::map<std::string, std::shared_ptr<VTKHandle>> m_vtk_auto;
  /** VTK writers that are executed manually */
  std::map<std::string, std::shared_ptr<VTKHandle>> m_vtk_manual;

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

  std::unique_ptr<pystencils::ContinuityKernel> m_continuity;

  // ResetFlux + external force
  // TODO: kernel for that
  // std::shared_ptr<ResetForce<PdfField, VectorField>> m_reset_force;

  [[nodiscard]] size_t stencil_size() const override { return FluxCount; }

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
  EKinWalberlaImpl(std::shared_ptr<LatticeWalberla> lattice,
                   FloatType diffusion, FloatType kT, FloatType valency,
                   Utils::Vector<FloatType, 3> ext_efield, FloatType density,
                   bool advection, bool friction_coupling)
      : EKinWalberlaBase<FloatType>(diffusion, kT, valency, ext_efield,
                                    advection, friction_coupling),
        m_lattice(std::move(lattice)) {
    m_density_field_id = field::addToStorage<DensityField>(
        m_lattice->get_blocks(), "density field", density, field::fzyx,
        m_lattice->get_ghost_layers());
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

    m_continuity = std::make_unique<pystencils::ContinuityKernel>(
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
  };

  void ghost_communication() override { (*m_full_communication)(); };

private:
  inline void kernel_boundary_density() {
    for (auto &block : *m_lattice->get_blocks()) {
      (*m_boundary_density)(&block);
    }
  }

  inline void kernel_boundary_flux() {
    for (auto &block : *m_lattice->get_blocks()) {
      (*m_boundary_flux)(&block);
    }
  }

  inline void kernel_continuity() {
    for (auto &block : *m_lattice->get_blocks()) {
      (*m_continuity).run(&block);
    }
  }

  inline void kernel_diffusion() {
    auto kernel = pystencils::DiffusiveFluxKernel(get_diffusion(),
                                                  m_flux_field_flattened_id,
                                                  m_density_field_flattened_id);

    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  inline void kernel_advection(const std::size_t &velocity_id) {
    auto kernel = pystencils::AdvectiveFluxKernel(m_flux_field_flattened_id,
                                                  m_density_field_id,
                                                  BlockDataID(velocity_id));
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  inline void kernel_friction_coupling(const std::size_t &force_id) {
    auto kernel = pystencils::FrictionCouplingKernel(
        get_diffusion(), BlockDataID(force_id), m_flux_field_flattened_id,
        get_kT());
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  inline void kernel_diffusion_electrostatic(const std::size_t &potential_id) {
    const auto ext_field = get_ext_efield();
    auto kernel = pystencils::DiffusiveFluxKernelWithElectrostatic(
        get_diffusion(), m_flux_field_flattened_id, BlockDataID(potential_id),
        m_density_field_flattened_id, ext_field[0], ext_field[1], ext_field[2],
        get_kT(), get_valency());
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  inline void kernel_migration() {}

public:
  void integrate(const std::size_t &potential_id,
                 const std::size_t &velocity_id,
                 const std::size_t &force_id) override {
    // early-breakout, has to be removed when reactions are included
    if (get_diffusion() == 0.)
      return;

    kernel_boundary_density();

    if (get_valency() != 0.) {
      if (potential_id == BlockDataID{}) {
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
      if (force_id == BlockDataID{}) {
        throw std::runtime_error("Walberla EK: friction coupling enabled but "
                                 "no force field accessible. force_id is " +
                                 std::to_string(force_id) +
                                 ". Have you enabled LB?");
      }
      kernel_friction_coupling(force_id);
    }

    if (get_advection()) {
      if (velocity_id == BlockDataID{}) {
        throw std::runtime_error("Walberla EK: advection enabled but no "
                                 "velocity field accessible. velocity_id is " +
                                 std::to_string(velocity_id) +
                                 ". Have you enabled LB?");
      }
      kernel_advection(velocity_id);
    }
    kernel_continuity();
    ghost_communication();

    // Handle VTK writers
    integrate_vtk_writers();
  };

  inline void integrate_vtk_writers() {
    for (const auto &it : m_vtk_auto) {
      auto &vtk_handle = it.second;
      if (vtk_handle->enabled) {
        vtk::writeFiles(vtk_handle->ptr)();
        vtk_handle->execution_count++;
      }
    }
  }

  [[nodiscard]] walberla::BlockDataID get_density_id() const override {
    return m_density_field_id;
  }

  // Density
  bool set_node_density(const Utils::Vector3i &node,
                        FloatType density) override {
    auto bc = get_block_and_cell(get_lattice(), node, false);
    if (!bc)
      return false;

    auto density_field =
        (*bc).block->template getData<DensityField>(m_density_field_id);
    density_field->get((*bc).cell) = density;

    return true;
  };

  [[nodiscard]] boost::optional<FloatType>
  get_node_density(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(get_lattice(), node, false);

    if (!bc)
      return {boost::none};

    auto density_field =
        (*bc).block->template getData<DensityField>(m_density_field_id);

    return {density_field->get((*bc).cell)};
  };

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

    reallocate_flux_boundary_field();

    return true;
  };

  bool remove_node_from_flux_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_flux->remove_node_from_boundary(node, *bc);

    reallocate_flux_boundary_field();

    return true;
  }

  bool set_node_density_boundary(const Utils::Vector3i &node,
                                 double density) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_density->set_node_value_at_boundary(node, density, *bc);

    reallocate_density_boundary_field();

    return true;
  }

  bool remove_node_from_density_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(get_lattice(), node, true);
    if (!bc)
      return false;

    m_boundary_density->remove_node_from_boundary(node, *bc);

    reallocate_density_boundary_field();

    return true;
  }

  boost::optional<bool>
  get_node_is_flux_boundary(const Utils::Vector3i &node,
                            bool consider_ghosts) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    return {m_boundary_flux->node_is_boundary(*bc)};
  }

  boost::optional<bool>
  get_node_is_density_boundary(const Utils::Vector3i &node,
                               bool consider_ghosts) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    return {m_boundary_density->node_is_boundary(*bc)};
  }

  boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
    if (!bc)
      return {boost::none};

    return {m_boundary_density->node_is_boundary(*bc) or
            m_boundary_flux->node_is_boundary(*bc)};
  };

  void update_flux_boundary_from_shape(
      const std::vector<int> &raster_flat,
      const std::vector<double> &flux_flat) override {
    // reshape grids
    auto const grid_size = get_lattice().get_grid_dimensions();
    std::vector<Utils::Vector3d> flux_vectors =
        fill_3D_vector_array(grid_size, flux_flat);
    boost::const_multi_array_ref<Utils::Vector3d, 3> fluxes(flux_vectors.data(),
                                                            grid_size);
    boost::const_multi_array_ref<int, 3> raster(raster_flat.data(), grid_size);

    auto const &blocks = get_lattice().get_blocks();
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
      // lattice constant is 1
      auto const left = block->getAABB().min();
      auto const off_i = int_c(left[0]);
      auto const off_j = int_c(left[1]);
      auto const off_k = int_c(left[2]);

      // Get field data which knows about the indices
      // In the loop, x,y,z are in block-local coordinates
      auto const n_ghost_layers = get_lattice().get_ghost_layers();
      auto const ghosts = static_cast<int>(n_ghost_layers);
      auto density_field =
          block->template getData<DensityField>(m_density_field_id);
      for (int i = -ghosts; i < int_c(density_field->xSize() + ghosts); ++i) {
        for (int j = -ghosts; j < int_c(density_field->ySize() + ghosts); ++j) {
          for (int k = -ghosts; k < int_c(density_field->zSize() + ghosts);
               ++k) {
            Utils::Vector3i const node{{off_i + i, off_j + j, off_k + k}};
            auto const idx = (node + grid_size) % grid_size;
            if (raster(idx)) {
              auto const bc = get_block_and_cell(get_lattice(), node, true);
              auto const &vel = fluxes(idx);
              m_boundary_flux->set_node_value_at_boundary(node, vel, *bc);
            }
          }
        }
      }
    }
    reallocate_flux_boundary_field();
  }

  void
  update_flux_boundary_from_list(std::vector<int> const &nodes_flat,
                                 std::vector<double> const &vel_flat) override {
    assert(nodes_flat.size() == vel_flat.size());
    assert(nodes_flat.size() % 3u == 0);
    for (std::size_t i = 0; i < nodes_flat.size(); i += 3) {
      auto const node = Utils::Vector3i(&nodes_flat[i], &nodes_flat[i + 3]);
      auto const vel = Utils::Vector3d(&vel_flat[i], &vel_flat[i + 3]);
      auto const bc = get_block_and_cell(get_lattice(), node, true);
      if (bc) {
        m_boundary_flux->set_node_value_at_boundary(node, vel, *bc);
      }
    }
    reallocate_flux_boundary_field();
  }

  void update_density_boundary_from_shape(
      const std::vector<int> &raster_flat,
      const std::vector<double> &density_flat) override {
    // reshape grids
    auto const grid_size = get_lattice().get_grid_dimensions();
    std::vector<double> density_vector =
        fill_3D_scalar_array(grid_size, density_flat);
    boost::const_multi_array_ref<double, 3> densities(density_vector.data(),
                                                      grid_size);
    boost::const_multi_array_ref<int, 3> raster(raster_flat.data(), grid_size);

    auto const &blocks = get_lattice().get_blocks();
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {
      // lattice constant is 1
      auto const left = block->getAABB().min();
      auto const off_i = int_c(left[0]);
      auto const off_j = int_c(left[1]);
      auto const off_k = int_c(left[2]);

      // Get field data which knows about the indices
      // In the loop, x,y,z are in block-local coordinates
      auto const n_ghost_layers = get_lattice().get_ghost_layers();
      auto const ghosts = static_cast<int>(n_ghost_layers);
      auto density_field =
          block->template getData<DensityField>(m_density_field_id);
      for (int i = -ghosts; i < int_c(density_field->xSize() + ghosts); ++i) {
        for (int j = -ghosts; j < int_c(density_field->ySize() + ghosts); ++j) {
          for (int k = -ghosts; k < int_c(density_field->zSize() + ghosts);
               ++k) {
            Utils::Vector3i const node{{off_i + i, off_j + j, off_k + k}};
            auto const idx = (node + grid_size) % grid_size;
            if (raster(idx)) {
              auto const bc = get_block_and_cell(get_lattice(), node, true);
              auto const &val = densities(idx);
              m_boundary_density->set_node_value_at_boundary(node, val, *bc);
            }
          }
        }
      }
    }
    reallocate_density_boundary_field();
  }

  void update_density_boundary_from_list(
      std::vector<int> const &nodes_flat,
      std::vector<double> const &density_flat) override {
    assert(nodes_flat.size() == vel_flat.size());
    assert(nodes_flat.size() % 3u == 0);
    for (std::size_t i = 0; i < nodes_flat.size(); i += 1) {
      auto const node =
          Utils::Vector3i(&nodes_flat[3 * i], &nodes_flat[3 * i + 3]);
      auto const val = density_flat[i];
      auto const bc = get_block_and_cell(get_lattice(), node, true);
      if (bc) {
        m_boundary_density->set_node_value_at_boundary(node, val, *bc);
      }
    }
    reallocate_density_boundary_field();
  }

  void reallocate_flux_boundary_field() { m_boundary_flux->boundary_update(); }

  void reallocate_density_boundary_field() {
    m_boundary_density->boundary_update();
  }

  [[nodiscard]] uint64_t get_rng_state() const override {
    throw std::runtime_error("The LB does not use a random number generator");
  };
  void set_rng_state(uint64_t counter) override {
    throw std::runtime_error("The LB does not use a random number generator");
  };

  [[nodiscard]] LatticeWalberla &get_lattice() const override {
    return *m_lattice;
  };

  // Grid, domain, halo

  [[nodiscard]] std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>>
  node_indices_positions(bool include_ghosts = false) const override {
    int ghost_offset = 0;
    if (include_ghosts)
      ghost_offset = m_lattice->get_ghost_layers();
    std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>> res;
    for (auto block = m_lattice->get_blocks()->begin();
         block != m_lattice->get_blocks()->end(); ++block) {
      auto left = block->getAABB().min();
      // Lattice constant is 1, node centers are offset by .5
      Utils::Vector3d pos_offset =
          to_vector3d(left) + Utils::Vector3d::broadcast(.5);

      // Lattice constant is 1, so cast left corner position to ints
      Utils::Vector3i index_offset =
          Utils::Vector3i{int(left[0]), int(left[1]), int(left[2])};

      // Get field data which knows about the indices
      // In the loop, x,y,z are in block-local coordinates
      auto density_field =
          block->template getData<DensityField>(m_density_field_id);
      for (int x = -ghost_offset;
           x < int(density_field->xSize()) + ghost_offset; x++) {
        for (int y = -ghost_offset;
             y < int(density_field->ySize()) + ghost_offset; y++) {
          for (int z = -ghost_offset;
               z < int(density_field->zSize()) + ghost_offset; z++) {
            res.emplace_back(
                index_offset + Utils::Vector3i{x, y, z},
                pos_offset + Utils::Vector3d{double(x), double(y), double(z)});
          }
        }
      }
    }
    return res;
  };

  class vtk_runtime_error : public std::runtime_error {
  public:
    explicit vtk_runtime_error(std::string const &vtk_uid,
                               std::string const &reason)
        : std::runtime_error("EKVTKOutput object '" + vtk_uid + "' " + reason) {
    }
  };

  [[nodiscard]] std::shared_ptr<VTKHandle>
  create_vtk(int delta_N, int initial_count, int flag_observables,
             std::string const &identifier, std::string const &base_folder,
             std::string const &prefix) override {
    // VTKOutput object must be unique
    std::stringstream unique_identifier;
    unique_identifier << base_folder << "/" << identifier;
    std::string const vtk_uid = unique_identifier.str();
    if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end() or
        m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
      throw vtk_runtime_error(vtk_uid, "already exists");
    }

    // instantiate VTKOutput object
    auto const &blocks = get_lattice().get_blocks();
    auto const write_freq = (delta_N) ? static_cast<unsigned int>(delta_N) : 1u;
    auto density_field_vtk = vtk::createVTKOutput_BlockData(
        blocks, identifier, uint_c(write_freq), uint_c(0), false, base_folder,
        prefix, true, true, true, true, uint_c(initial_count));
    field::FlagFieldCellFilter<FlagField> fluid_filter(m_flag_field_density_id);
    fluid_filter.addFlag(Boundary_flag);
    density_field_vtk->addCellExclusionFilter(fluid_filter);

    // add writers
    if (flag_observables & static_cast<int>(EKOutputVTK::density)) {
      density_field_vtk->addCellDataWriter(
          make_shared<field::VTKWriter<DensityField, float>>(m_density_field_id,
                                                             "density"));
    }

    // register object
    auto vtk_handle =
        std::make_shared<VTKHandle>(density_field_vtk, initial_count, true);
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

  ~EKinWalberlaImpl() override = default;
};
} // namespace walberla

#endif // ESPRESSO_EKINWALBERLAIMPL_HPP
