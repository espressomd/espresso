#ifndef ESPRESSO_EKINWALBERLAIMPL_HPP
#define ESPRESSO_EKINWALBERLAIMPL_HPP

#include "blockforest/communication/UniformBufferedScheme.h"
#include "boundary/BoundaryHandling.h"
#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/lattice_model/D3Q27.h"
#include "timeloop/SweepTimeloop.h"

#include "EKinWalberlaBase.hpp"
#include "LatticeWalberla.hpp"
#include "walberla_utils.hpp"

#include "generated_kernels/AdvectiveFluxKernel.h"
#include "generated_kernels/ContinuityKernel.h"
#include "generated_kernels/DiffusiveFluxKernel.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostatic.h"
#include "generated_kernels/FrictionCouplingKernel.h"
#include "generated_kernels/NoFlux.h"

#include <memory>

namespace walberla {

// Flags marking fluid and boundaries
const FlagUID domain_flag("domain");
const FlagUID noflux_flag("noflux");

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

  /** VTK writers that are executed automatically */
  std::map<std::string, std::pair<std::shared_ptr<vtk::VTKOutput>, bool>>
      m_vtk_auto;
  /** VTK writers that are executed manually */
  std::map<std::string, std::shared_ptr<vtk::VTKOutput>> m_vtk_manual;

  // Block data access handles
  BlockDataID m_density_field_id;
  BlockDataID m_density_field_flattened_id;
  BlockDataID m_flag_field_id;

  BlockDataID m_flux_field_id;
  BlockDataID m_flux_field_flattened_id;

  /** Block forest */
  std::shared_ptr<LatticeWalberla> m_lattice;

  std::unique_ptr<pystencils::NoFlux> m_noflux;
  bool m_noflux_dirty = false;

  std::unique_ptr<pystencils::ContinuityKernel> m_continuity;

  // ResetFlux + external force
  // TODO: kernel for that
  // std::shared_ptr<ResetForce<PdfField, VectorField>> m_reset_force;

  [[nodiscard]] size_t stencil_size() const override { return FluxCount; }

  // Boundary handling
  // TODO: Boundary Handling for density and fluxes
  //  class BoundaryHandling {
  //  public:
  //    BoundaryHandling(const BlockDataID &flag_field_id,
  //                     const BlockDataID &flux_field_id)
  //        : m_flag_field_id(flag_field_id), m_flux_field_id(flux_field_id) {}
  //
  //    Boundaries *operator()(IBlock *const block) {
  //
  //      auto *flag_field = block->template
  //      getData<FlagField>(m_flag_field_id); auto *pdf_field = block->template
  //      getData<FluxField>(m_flux_field_id);
  //
  //      const auto domain = flag_field->flagExists(domain_flag)
  //                              ? flag_field->getFlag(domain_flag)
  //                              : flag_field->registerFlag(domain_flag);
  //
  //      pystencils::NoFlux noflux(block, m_flux_field_id);
  //
  //      return new Boundaries(
  //          "boundary handling", flag_field, fluid,
  //          UBB("velocity bounce back", UBB_flag, pdf_field, nullptr));
  //    }
  //
  //  private:
  //    const BlockDataID m_flag_field_id;
  //    const BlockDataID m_flux_field_id;
  //  };

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

    m_noflux = std::make_unique<pystencils::NoFlux>(m_lattice->get_blocks(),
                                                    m_flux_field_flattened_id);

    m_continuity = std::make_unique<pystencils::ContinuityKernel>(
        m_flux_field_flattened_id, m_density_field_flattened_id);

    // Init and register flag field (domain/boundary)
    m_flag_field_id = field::addFlagFieldToStorage<FlagField>(
        m_lattice->get_blocks(), "flag field", m_lattice->get_ghost_layers());

    for (auto block = m_lattice->get_blocks()->begin();
         block != m_lattice->get_blocks()->end(); ++block) {
      auto *flagField = block->template getData<FlagField>(m_flag_field_id);

      flagField->registerFlag(domain_flag);
      flagField->registerFlag(noflux_flag);
    }

    // set domain flag everywhere
    clear_boundaries();

    // m_noflux = pystencils::NoFlux(m_lattice,
    // m_flux_field_flattened_id);

    // TODO: boundary handlings
    // Register boundary handling
    //    m_boundary_handling_id =
    //        m_blockforest->get_blocks()->addBlockData<Boundaries>(
    //            LBBoundaryHandling(m_flag_field_id, m_pdf_field_id),
    //            "boundary handling");
    //    clear_boundaries();

    m_full_communication =
        std::make_shared<FullCommunicator>(m_lattice->get_blocks());
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<DensityField>>(
            m_density_field_id));

    // reset flux field
    // diffusive flux + migrative flux + output friction force if field is
    // provided... add advective flux continuum update

    //    m_time_loop->add() << timeloop::Sweep(makeSharedSweep(m_reset_force),
    //                                          "Reset force fields");
    //    m_time_loop->add() << timeloop::Sweep(collide, "LB collide")
    //                       << timeloop::AfterFunction(
    //                           *m_pdf_streaming_communication,
    //                           "communication");
    //    m_time_loop->add() << timeloop::Sweep(
    //        Boundaries::getBlockSweep(m_boundary_handling_id), "boundary
    //        handling");
    //    m_time_loop->add() << timeloop::Sweep(stream, "LB stream")
    //                       << timeloop::AfterFunction(*m_full_communication,
    //                                                  "communication");

    // Synchronize ghost layers
    (*m_full_communication)();
  };

  void ghost_communication() override { (*m_full_communication)(); };

private:
  inline void kernel_noflux() {
    for (auto &block : *m_lattice->get_blocks()) {
      (*m_noflux).run(&block);
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

  inline void kernel_advection(const BlockDataID &velocity_id) {
    auto kernel = pystencils::AdvectiveFluxKernel(
        m_flux_field_flattened_id, m_density_field_id, velocity_id);
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  inline void kernel_friction_coupling(const BlockDataID &force_id) {
    auto kernel = pystencils::FrictionCouplingKernel(
        get_diffusion(), force_id, m_flux_field_flattened_id, get_kT());
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  inline void kernel_diffusion_electrostatic(const BlockDataID &potential_id) {
    const auto ext_field = get_ext_efield();
    auto kernel = pystencils::DiffusiveFluxKernelWithElectrostatic(
        get_diffusion(), m_flux_field_flattened_id, potential_id,
        m_density_field_flattened_id, ext_field[0], ext_field[1], ext_field[2],
        get_kT(), get_valency());
    for (auto &block : *m_lattice->get_blocks()) {
      kernel.run(&block);
    }
  }

  inline void kernel_migration() {}

public:
  void integrate(const BlockDataID &potential_id,
                 const BlockDataID &velocity_id,
                 const BlockDataID &force_id) override {
    if (m_noflux_dirty) {
      boundary_noflux_update();
      m_noflux_dirty = false;
    }

    // early-breakout, has to be removed when reactions are included
    if (get_diffusion() == 0.)
      return;

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
    kernel_noflux();
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
    for (const auto &it : m_vtk_auto) {
      if (it.second.second)
        vtk::writeFiles(it.second.first)();
    }
  };

  [[nodiscard]] walberla::BlockDataID get_density_id() const override {
    return m_density_field_id;
  }

  // Density
  bool set_node_density(const Utils::Vector3i &node,
                        FloatType density) override {
    auto bc = get_block_and_cell(node, false, m_lattice->get_blocks(),
                                 m_lattice->get_ghost_layers());
    if (!bc)
      return false;

    auto density_field =
        (*bc).block->template getData<DensityField>(m_density_field_id);
    density_field->get((*bc).cell) = density;

    return true;
  };

  [[nodiscard]] boost::optional<FloatType>
  get_node_density(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, false, m_lattice->get_blocks(),
                                 m_lattice->get_ghost_layers());

    if (!bc)
      return {boost::none};

    auto density_field =
        (*bc).block->template getData<DensityField>(m_density_field_id);

    return {density_field->get((*bc).cell)};
  };

  void boundary_noflux_update() {
    m_noflux.get()->template fillFromFlagField<FlagField>(
        m_lattice->get_blocks(), m_flag_field_id, noflux_flag, domain_flag);
  }

  [[nodiscard]] bool
  set_node_noflux_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(node, true, m_lattice->get_blocks(),
                                 m_lattice->get_ghost_layers());
    if (!bc)
      return false;

    auto *flagField = (*bc).block->template getData<FlagField>(m_flag_field_id);
    flagField->addFlag((*bc).cell, flagField->getFlag(noflux_flag));

    m_noflux_dirty = true;
    return true;
  }

  [[nodiscard]] bool
  remove_node_from_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(node, true, m_lattice->get_blocks(),
                                 m_lattice->get_ghost_layers());
    if (!bc)
      return false;

    auto *flagField = (*bc).block->template getData<FlagField>(m_flag_field_id);
    flagField->removeFlag((*bc).cell, flagField->getFlag(noflux_flag));

    m_noflux_dirty = true;
    return true;
  };

  [[nodiscard]] boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(node, consider_ghosts, m_lattice->get_blocks(),
                                 m_lattice->get_ghost_layers());
    if (!bc)
      return {boost::none};

    auto *flagField = (*bc).block->template getData<FlagField>(m_flag_field_id);
    // TODO: this has to be adapted when other BCs are added
    return {flagField->isPartOfMaskSet((*bc).cell,
                                       flagField->getFlag(noflux_flag))};
  };

  void clear_boundaries() override {
    const CellInterval &domain_bb_in_global_cell_coordinates =
        m_lattice->get_blocks()->getCellBBFromAABB(
            m_lattice->get_blocks()->begin()->getAABB().getExtended(
                FloatType(m_lattice->get_ghost_layers())));
    for (auto block = m_lattice->get_blocks()->begin();
         block != m_lattice->get_blocks()->end(); ++block) {

      auto *flagField = block->template getData<FlagField>(m_flag_field_id);

      CellInterval domain_bb(domain_bb_in_global_cell_coordinates);
      m_lattice->get_blocks()->transformGlobalToBlockLocalCellInterval(
          domain_bb, *block);

      for (auto cell : domain_bb) {
        flagField->addFlag(cell, flagField->getFlag(domain_flag));
      }
    }

    m_noflux_dirty = true;
  };

  [[nodiscard]] uint64_t get_rng_state() const override {
    throw std::runtime_error("The LB does not use a random number generator");
  };
  void set_rng_state(uint64_t counter) override {
    throw std::runtime_error("The LB does not use a random number generator");
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

  void create_vtk(unsigned delta_N, unsigned initial_count,
                  unsigned flag_observables, std::string const &identifier,
                  std::string const &base_folder,
                  std::string const &prefix) override {
    // VTKOuput object must be unique
    std::stringstream unique_identifier;
    unique_identifier << base_folder << "/" << identifier;
    std::string const vtk_uid = unique_identifier.str();
    if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end() or
        m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " already exists");
    }

    // instantiate VTKOutput object
    unsigned const write_freq = (delta_N) ? static_cast<unsigned>(delta_N) : 1u;
    auto vtk_writer = vtk::createVTKOutput_BlockData(
        m_lattice->get_blocks(), identifier, uint_c(write_freq), uint_c(0),
        false, base_folder, prefix, true, true, true, true,
        uint_c(initial_count));
    field::FlagFieldCellFilter<FlagField> domain_filter(m_flag_field_id);
    domain_filter.addFlag(domain_flag);
    vtk_writer->addCellInclusionFilter(domain_filter);

    // add writers
    if (static_cast<unsigned>(EKOutputVTK::density) & flag_observables) {
      vtk_writer->addCellDataWriter(
          make_shared<field::VTKWriter<DensityField, float>>(m_density_field_id,
                                                             "EKDensity"));
    }

    // register object
    if (delta_N) {
      m_vtk_auto[vtk_uid] = {vtk_writer, true};
    } else {
      m_vtk_manual[vtk_uid] = vtk_writer;
    }
  }

  /** Manually call a VTK callback */
  void write_vtk(std::string const &vtk_uid) override {
    if (m_vtk_auto.find(vtk_uid) != m_vtk_auto.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " is an automatic observable");
    }
    if (m_vtk_manual.find(vtk_uid) == m_vtk_manual.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " doesn't exist");
    }
    vtk::writeFiles(m_vtk_manual[vtk_uid])();
  }

  /** Activate or deactivate a VTK callback */
  void switch_vtk(std::string const &vtk_uid, int status) override {
    if (m_vtk_manual.find(vtk_uid) != m_vtk_manual.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " is a manual observable");
    }
    if (m_vtk_auto.find(vtk_uid) == m_vtk_auto.end()) {
      throw std::runtime_error("VTKOutput object " + vtk_uid +
                               " doesn't exist");
    }
    m_vtk_auto[vtk_uid].second = status;
  }

  ~EKinWalberlaImpl() override = default;
};
} // namespace walberla

#endif // ESPRESSO_EKINWALBERLAIMPL_HPP
