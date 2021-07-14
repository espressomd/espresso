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
#include "walberla_utils.hpp"

#include "generated_kernels/ContinuityKernel.h"
#include "generated_kernels/DiffusiveFluxKernel.h"
#include "generated_kernels/NoFlux.h"

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
  const WalberlaBlockForest *m_blockforest;

  pystencils::NoFlux *m_noflux;

  std::shared_ptr<timeloop::SweepTimeloop> m_time_loop;

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
  EKinWalberlaImpl(const WalberlaBlockForest *blockforest, FloatType diffusion,
                   FloatType kT, FloatType density)
      : EKinWalberlaBase<FloatType>(diffusion, kT), m_blockforest{blockforest} {
    m_density_field_id = field::addToStorage<DensityField>(
        get_blockforest()->get_blocks(), "density field", density, field::fzyx,
        get_blockforest()->get_ghost_layers());
    m_density_field_flattened_id =
        field::addFlattenedShallowCopyToStorage<DensityField>(
            get_blockforest()->get_blocks(), m_density_field_id,
            "flattened density field");
    m_flux_field_id = field::addToStorage<FluxField>(
        get_blockforest()->get_blocks(), "flux field", FloatType{0},
        field::fzyx, get_blockforest()->get_ghost_layers());
    m_flux_field_flattened_id =
        field::addFlattenedShallowCopyToStorage<FluxField>(
            get_blockforest()->get_blocks(), m_flux_field_id,
            "flattened flux field");

    m_noflux = new pystencils::NoFlux(get_blockforest()->get_blocks(),
                                      m_flux_field_flattened_id);

    // Init and register flag field (domain/boundary)
    m_flag_field_id = field::addFlagFieldToStorage<FlagField>(
        get_blockforest()->get_blocks(), "flag field",
        get_blockforest()->get_ghost_layers());

    for (auto block = get_blockforest()->get_blocks()->begin();
         block != get_blockforest()->get_blocks()->end(); ++block) {
      auto *flagField = block->template getData<FlagField>(m_flag_field_id);

      flagField->registerFlag(domain_flag);
      flagField->registerFlag(noflux_flag);
    }

    // set domain flag everywhere
    clear_boundaries();

    // m_noflux = pystencils::NoFlux(get_blockforest(),
    // m_flux_field_flattened_id);

    // TODO: boundary handlings
    // Register boundary handling
    //    m_boundary_handling_id =
    //        m_blockforest->get_blocks()->addBlockData<Boundaries>(
    //            LBBoundaryHandling(m_flag_field_id, m_pdf_field_id),
    //            "boundary handling");
    //    clear_boundaries();

    m_full_communication =
        std::make_shared<FullCommunicator>(get_blockforest()->get_blocks());
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<DensityField>>(
            m_density_field_id));

    // TODO: kernels references

    // Add steps to the integration loop
    m_time_loop = std::make_shared<timeloop::SweepTimeloop>(
        get_blockforest()->get_blocks()->getBlockStorage(), 1);
    // TODO: figure out if there is a less-hacky solution to this because this
    // re-instantiates the Diffusive-Flux kernel for every call...
    m_time_loop->add() << Sweep(
        [this](IBlock *block) {
          return pystencils::DiffusiveFluxKernel(this->get_diffusion(),
                                                 m_flux_field_flattened_id,
                                                 m_density_field_flattened_id)
              .run(block);
        },
        "ekin diffusive flux");
    m_time_loop->add() << Sweep(*m_noflux, "no flux boundary condition");
    m_time_loop->add()
        << Sweep(pystencils::ContinuityKernel(m_flux_field_flattened_id,
                                              m_density_field_flattened_id),
                 "ekin continuity")
        << AfterFunction(*m_full_communication, "ekin density communication");

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

  void integrate() override {
    m_time_loop->singleStep();

    // Handle VTK writers
    for (const auto &it : m_vtk_auto) {
      if (it.second.second)
        vtk::writeFiles(it.second.first)();
    }
  };

  [[nodiscard]] const WalberlaBlockForest *get_blockforest() const override {
    return m_blockforest;
  };

  // Density
  bool set_node_density(const Utils::Vector3i &node,
                        FloatType density) override {
    auto bc = get_block_and_cell(node, false, get_blockforest()->get_blocks(),
                                 get_blockforest()->get_ghost_layers());
    if (!bc)
      return false;

    auto density_field =
        (*bc).block->template getData<DensityField>(m_density_field_id);
    density_field->get((*bc).cell) = density;

    return true;
  };

  [[nodiscard]] boost::optional<FloatType>
  get_node_density(const Utils::Vector3i &node) const override {
    auto bc = get_block_and_cell(node, false, get_blockforest()->get_blocks(),
                                 get_blockforest()->get_ghost_layers());
    if (!bc)
      return {boost::none};

    auto density_field =
        (*bc).block->template getData<DensityField>(m_density_field_id);

    return {density_field->get((*bc).cell)};
  };

  void boundary_noflux_update() {
    m_noflux->template fillFromFlagField<FlagField>(
        get_blockforest()->get_blocks(), m_flag_field_id, noflux_flag,
        domain_flag);
  }

  [[nodicard]] bool
  set_node_noflux_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(node, true, get_blockforest()->get_blocks(),
                                 get_blockforest()->get_ghost_layers());
    if (!bc)
      return false;

    auto *flagField = (*bc).block->template getData<FlagField>(m_flag_field_id);
    flagField->addFlag((*bc).cell, flagField->getFlag(noflux_flag));

    boundary_noflux_update();
    return true;
  }

  [[nodiscard]] bool
  remove_node_from_boundary(const Utils::Vector3i &node) override {
    auto bc = get_block_and_cell(node, true, get_blockforest()->get_blocks(),
                                 get_blockforest()->get_ghost_layers());
    if (!bc)
      return false;

    auto *flagField = (*bc).block->template getData<FlagField>(m_flag_field_id);
    flagField->removeFlag((*bc).cell, flagField->getFlag(noflux_flag));

    boundary_noflux_update();
    return true;
  };

  [[nodiscard]] boost::optional<bool>
  get_node_is_boundary(const Utils::Vector3i &node,
                       bool consider_ghosts = false) const override {
    auto bc = get_block_and_cell(node, consider_ghosts,
                                 get_blockforest()->get_blocks(),
                                 get_blockforest()->get_ghost_layers());
    if (!bc)
      return {boost::none};

    auto *flagField = (*bc).block->template getData<FlagField>(m_flag_field_id);
    return {!flagField->isFlagSet((*bc).cell, flagField->getFlag(domain_flag))};
  };

  void clear_boundaries() override {
    const CellInterval &domain_bb_in_global_cell_coordinates =
        get_blockforest()->get_blocks()->getCellBBFromAABB(
            get_blockforest()->get_blocks()->begin()->getAABB().getExtended(
                FloatType(get_blockforest()->get_ghost_layers())));
    for (auto block = get_blockforest()->get_blocks()->begin();
         block != get_blockforest()->get_blocks()->end(); ++block) {

      auto *flagField = block->template getData<FlagField>(m_flag_field_id);

      CellInterval domain_bb(domain_bb_in_global_cell_coordinates);
      get_blockforest()->get_blocks()->transformGlobalToBlockLocalCellInterval(
          domain_bb, *block);

      for (auto cell : domain_bb) {
        flagField->addFlag(cell, flagField->getFlag(domain_flag));
      }
    }

    boundary_noflux_update();
  };

  [[nodiscard]] uint64_t get_rng_state() const override {
    throw std::runtime_error("The LB does not use a random number generator");
  };
  void set_rng_state(uint64_t counter) override {
    throw std::runtime_error("The LB does not use a random number generator");
  };

  // Grid, domain, halo

  std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>>
  node_indices_positions(bool include_ghosts = false) const override {
    int ghost_offset = 0;
    if (include_ghosts)
      ghost_offset = get_blockforest()->get_ghost_layers();
    std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>> res;
    for (auto block = get_blockforest()->get_blocks()->begin();
         block != get_blockforest()->get_blocks()->end(); ++block) {
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
        get_blockforest()->get_blocks(), identifier, uint_c(write_freq),
        uint_c(0), false, base_folder, prefix, true, true, true, true,
        uint_c(initial_count));
    field::FlagFieldCellFilter<FlagField> domain_filter(m_flag_field_id);
    domain_filter.addFlag(domain_flag);
    vtk_writer->addCellInclusionFilter(domain_filter);

    // add writers
    if (static_cast<unsigned>(EKOutputVTK::density) & flag_observables) {
      // TODO: density writer...
      //      vtk_writer->addCellDataWriter(
      //          make_shared<lbm::DensityVTKWriter<LatticeModel, float>>(
      //              m_pdf_field_id, "DensityFromPDF"));
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
