#ifndef LB_WALBERLA_H
#define LB_WALBERLA_H

#include "config.hpp"

#ifdef LB_WALBERLA
#include "utils/Vector.hpp"

#undef PI
#include "blockforest/StructuredBlockForest.h"
#include "boost/optional.hpp"
#include "boost/tuple/tuple.hpp"
#include "boundary/BoundaryHandling.h"
#include "core/mpi/Environment.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/distributors/KernelDistributor.h"
#include "field/interpolators/TrilinearFieldInterpolator.h"
#include "lbm/boundary/UBB.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "timeloop/SweepTimeloop.h"

const walberla::FlagUID Fluid_flag("fluid");
const walberla::FlagUID UBB_flag("velocity bounce back");

inline Utils::Vector3d
to_vector3d(const walberla::Vector3<walberla::real_t> v) {
  return Utils::Vector3d{v[0], v[1], v[2]};
}
inline walberla::Vector3<walberla::real_t> to_vector3(const Utils::Vector3d v) {
  return walberla::Vector3<walberla::real_t>{v[0], v[1], v[2]};
}
inline Utils::Vector6d
to_vector6d(const walberla::Matrix3<walberla::real_t> m) {
  return Utils::Vector6d{m[0], m[3], m[4], m[6], m[7], m[8]};
}

template <typename PdfField_T, typename ForceField_T,
          typename BoundaryHandling_T>
class ResetForce {
public:
  ResetForce(walberla::BlockDataID pdf_field_id,
             walberla::BlockDataID force_field_id,
             walberla::BlockDataID force_field_from_md_id,
             walberla::BlockDataID boundary_handling_id)
      : m_pdf_field_id(pdf_field_id), m_force_field_id(force_field_id),
        m_force_field_from_md_id(force_field_from_md_id),
        m_boundary_handling_id(boundary_handling_id),
        m_ext_force(walberla::Vector3<walberla::real_t>{0, 0, 0}) {}

  void set_ext_force(const Utils::Vector3d &ext_force) {
    m_ext_force = to_vector3(ext_force);
  }

  Utils::Vector3d get_ext_force() const { return to_vector3d(m_ext_force); };

  void operator()(walberla::IBlock *block) {
    PdfField_T *pdf_field = block->getData<PdfField_T>(m_pdf_field_id);
    ForceField_T *force_field = block->getData<ForceField_T>(m_force_field_id);
    ForceField_T *force_field_from_md =
        block->getData<ForceField_T>(m_force_field_from_md_id);
    BoundaryHandling_T *boundary_handling =
        block->getData<BoundaryHandling_T>(m_boundary_handling_id);

    force_field->swapDataPointers(force_field_from_md);

    WALBERLA_FOR_ALL_CELLS_XYZ(force_field, {
      walberla::Cell cell(x, y, z);
      if (boundary_handling->isDomain(cell)) {
        force_field->get(cell) += m_ext_force;
        force_field->get(cell) /= pdf_field->getDensity(cell);
      }
    });
    WALBERLA_FOR_ALL_CELLS_XYZ(force_field_from_md, {
      walberla::Cell cell(x, y, z);
      if (boundary_handling->isDomain(cell)) {
        force_field_from_md->get(cell) = walberla::Vector3<walberla::real_t>{0};
      }
    });
  }

private:
  walberla::BlockDataID m_pdf_field_id, m_force_field_id,
      m_force_field_from_md_id, m_boundary_handling_id;
  walberla::Vector3<walberla::real_t> m_ext_force;
};
/** Class that runs and controls the LB on WaLBerla
 */
class LbWalberla {
  double m_skin;
  double m_agrid;
  walberla::uint_t n_ghost_layers() const {
    return walberla::uint_c(m_skin/m_agrid +1);
  }
  double m_tau;
  double m_density;           // initial density
  Utils::Vector3d m_velocity; // initial velocity
  Utils::Vector3i m_grid_dimensions;

  // Type definitions
  using vector_field_t =
      walberla::GhostLayerField<walberla::Vector3<walberla::real_t>, 1>;
  using force_model_t = walberla::lbm::force_model::GuoField<vector_field_t>;
  //  LatticeModel_T(lbm::collision_model::TRT::constructWithMagicNumber(
  //  walberla::lbm::collision_model::omegaFromViscosity(eta_s)),
  //  lbm::force_model::GuoField<ForceField_T>( forceFieldId ) );
  using Lattice_model_t =
      walberla::lbm::D3Q19<walberla::lbm::collision_model::TRT, false,
                           force_model_t>;
  using Flag_field_t = walberla::FlagField<walberla::uint8_t>;
  using Pdf_field_t = walberla::lbm::PdfField<Lattice_model_t>;

  using UBB_t = walberla::lbm::UBB<Lattice_model_t, walberla::uint8_t>;
  using Boundary_handling_t =
      walberla::BoundaryHandling<Flag_field_t, Lattice_model_t::Stencil, UBB_t>;

  std::shared_ptr<ResetForce<Pdf_field_t, vector_field_t, Boundary_handling_t>>
      m_reset_force;

  class LB_boundary_handling {
  public:
    LB_boundary_handling(const walberla::BlockDataID &flag_field_id,
                         const walberla::BlockDataID &pdf_field_id)
        : m_flag_field_id(flag_field_id), m_pdf_field_id(pdf_field_id) {}

    Boundary_handling_t *operator()(walberla::IBlock *const block) {

      Flag_field_t *flag_field = block->getData<Flag_field_t>(m_flag_field_id);
      Pdf_field_t *pdf_field = block->getData<Pdf_field_t>(m_pdf_field_id);

      const uint8_t fluid = flag_field->registerFlag(Fluid_flag);

      return new Boundary_handling_t(
          "boundary handling", flag_field, fluid,
          UBB_t("velocity bounce back", UBB_flag, pdf_field, nullptr));
    }

  private:
    const walberla::BlockDataID m_flag_field_id;
    const walberla::BlockDataID m_pdf_field_id;
  };

  struct BlockAndCell {
    walberla::IBlock *block;
    walberla::Cell cell;
  };

public:
  LbWalberla(double viscosity, double density, double agrid, double tau,
             const Utils::Vector3d &box_dimensions,
             const Utils::Vector3i &node_grid, double skin);

  void integrate();
  std::pair<Utils::Vector3d, Utils::Vector3d> get_local_domain() {
    // We only have one block per mpi rank
    assert(m_blocks->begin()++ == m_blocks->end());
    
    auto const ab=m_blocks->begin()->getAABB();
    return 
      {to_vector3d(ab.min()),to_vector3d(ab.max())};
  }
  boost::optional<Utils::Vector3d>
  get_node_velocity(const Utils::Vector3i node) const;
  bool set_node_velocity(const Utils::Vector3i &node, const Utils::Vector3d v);
  boost::optional<Utils::Vector19d>
  get_node_pop(const Utils::Vector3i node) const;
  bool set_node_pop(const Utils::Vector3i &node, const Utils::Vector19d pop);

  bool add_force_at_pos(const Utils::Vector3d &position,
                        const Utils::Vector3d &force);
  boost::optional<Utils::Vector3d>
  get_velocity_at_pos(const Utils::Vector3d &position) const;
  boost::optional<double> get_density_at_pos(const Utils::Vector3d &position);
  //  Utils::Vector3d get_stress_at_pos(const Utils::Vector3d& position);

  boost::optional<double> get_node_density(const Utils::Vector3i node) const;
  bool set_node_density(const Utils::Vector3i node, const double density);

  double get_density() const { return m_density; };

  boost::optional<Utils::Vector3d>
  get_node_velocity_at_boundary(const Utils::Vector3i &node) const;

  boost::optional<Utils::Vector6d>
  get_node_pressure_tensor(const Utils::Vector3i node) const;

  bool set_node_velocity_at_boundary(const Utils::Vector3i node,
                                     const Utils::Vector3d v);
  bool remove_node_from_boundary(const Utils::Vector3i &node);
  boost::optional<bool> get_node_is_boundary(const Utils::Vector3i &node) const;

  Utils::Vector3d get_momentum() const;

  void set_external_force(const Utils::Vector3d &ext_force) {
    m_reset_force->set_ext_force(ext_force);
  }
  Utils::Vector3d get_external_force() {
    return m_reset_force->get_ext_force();
  }

  void set_viscosity(double viscosity);
  double get_viscosity();

  Utils::Vector3i get_grid_dimensions() { return m_grid_dimensions; };

  double get_grid_spacing() { return m_agrid; }
  double get_tau() { return m_tau; }

  bool node_in_local_domain(const Utils::Vector3i &node) const;
  bool pos_in_local_domain(const Utils::Vector3d &pos) const;

  const Utils::Vector<walberla::uint_t, 19> es_pop_index_to_walberla_pop_index{
      Lattice_model_t::Stencil::idx[walberla::stencil::C],
      Lattice_model_t::Stencil::idx[walberla::stencil::E],
      Lattice_model_t::Stencil::idx[walberla::stencil::W],
      Lattice_model_t::Stencil::idx[walberla::stencil::N],
      Lattice_model_t::Stencil::idx[walberla::stencil::S],
      Lattice_model_t::Stencil::idx[walberla::stencil::T],
      Lattice_model_t::Stencil::idx[walberla::stencil::B],
      Lattice_model_t::Stencil::idx[walberla::stencil::NE],
      Lattice_model_t::Stencil::idx[walberla::stencil::SW],
      Lattice_model_t::Stencil::idx[walberla::stencil::SE],
      Lattice_model_t::Stencil::idx[walberla::stencil::NW],
      Lattice_model_t::Stencil::idx[walberla::stencil::TE],
      Lattice_model_t::Stencil::idx[walberla::stencil::BW],
      Lattice_model_t::Stencil::idx[walberla::stencil::BE],
      Lattice_model_t::Stencil::idx[walberla::stencil::TW],
      Lattice_model_t::Stencil::idx[walberla::stencil::TN],
      Lattice_model_t::Stencil::idx[walberla::stencil::BS],
      Lattice_model_t::Stencil::idx[walberla::stencil::BN],
      Lattice_model_t::Stencil::idx[walberla::stencil::TS]};

private:
  boost::optional<BlockAndCell>
  get_block_and_cell(const Utils::Vector3i &node, bool consider_ghost_layers=false) const;
  walberla::IBlock* get_block(const Utils::Vector3d &pos, bool consider_ghost_layers) const;
  walberla::BlockDataID m_pdf_field_id;
  walberla::BlockDataID m_flag_field_id;
  walberla::BlockDataID m_force_field_id;
  walberla::BlockDataID m_force_field_from_md_id;
  walberla::BlockDataID m_force_distributor_id;
  walberla::BlockDataID m_velocity_adaptor_id;
  walberla::BlockDataID m_velocity_interpolator_id;

  std::shared_ptr<walberla::blockforest::StructuredBlockForest> m_blocks;
  std::shared_ptr<walberla::timeloop::SweepTimeloop> m_time_loop;
  std::shared_ptr<walberla::mpi::Environment> m_env;
  walberla::BlockDataID m_boundary_handling_id;
  std::shared_ptr<Lattice_model_t> m_lattice_model;
  // Adaptors
  using DensityAdaptor = walberla::lbm::Adaptor<Lattice_model_t>::Density;
  using VelocityAdaptor =
      walberla::lbm::Adaptor<Lattice_model_t>::VelocityVector;

  using Vector_field_distributor_t =
      walberla::field::KernelDistributor<vector_field_t, Flag_field_t>;
  using VectorFieldAdaptorInterpolator =
      walberla::field::TrilinearFieldInterpolator<VelocityAdaptor,
                                                  Flag_field_t>;

  void empty_flags_for_boundary(
      std::shared_ptr<walberla::StructuredBlockForest> &blocks,
      const walberla::BlockDataID &boundary_handling_id) {
    using namespace walberla;
    const CellInterval &domain_bb_in_global_cell_coordinates =
        blocks->getDomainCellBB();
    for (auto block = blocks->begin(); block != blocks->end(); ++block) {

      Boundary_handling_t *boundary_handling =
          block->getData<Boundary_handling_t>(boundary_handling_id);

      CellInterval domain_bb(domain_bb_in_global_cell_coordinates);
      blocks->transformGlobalToBlockLocalCellInterval(domain_bb, *block);

      boundary_handling->fillWithDomain(domain_bb);
    }
  }
};

#endif // LB_WALBERLA

#endif // LB_WALBERLA_H
