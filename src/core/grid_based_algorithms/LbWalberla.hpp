#ifndef LB_WALBERLA_H
#define LB_WALBERLA_H

#ifdef LB_WALBERLA

#include "utils/Vector.hpp"
#undef PI
#include "blockforest/StructuredBlockForest.h"
#include "field/GhostLayerField.h"
#include "boundary/BoundaryHandling.h"
#include "core/mpi/Environment.h"
#include "field/FlagField.h"
#include "field/distributors/KernelDistributor.h"
#include "field/interpolators/TrilinearFieldInterpolator.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimpleUBB.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "timeloop/SweepTimeloop.h"
#include "vtk/VTKOutput.h"

const walberla::FlagUID Fluid_flag("fluid");
const walberla::FlagUID UBB_flag("velocity bounce back");
const walberla::FlagUID No_slip_flag("no slip");

/** Class that runs and controls the LB on WaLBerla
 */
class LbWalberla {
  double m_skin;
  double m_agrid;
  Vector3d m_ext_force;

  // Type definitions
  using vector_field_t =
      walberla::GhostLayerField<walberla::Vector3<walberla::real_t>, 1>;
  using force_model_t = walberla::lbm::force_model::GuoField<vector_field_t>;
  using Lattice_model_t =
      walberla::lbm::D3Q19<walberla::lbm::collision_model::SRT, false,
                           force_model_t>;
  using Flag_field_t = walberla::FlagField<walberla::uint8_t>;
  using Pdf_field_t = walberla::lbm::PdfField<Lattice_model_t>;

  using No_slip_t = walberla::lbm::NoSlip<Lattice_model_t, walberla::uint8_t>;
  using UBB_t = walberla::lbm::SimpleUBB<Lattice_model_t, walberla::uint8_t>;
  using Boundary_conditions_t = boost::tuples::tuple<No_slip_t, UBB_t>;
  using Boundary_handling_t =
      walberla::BoundaryHandling<Flag_field_t, Lattice_model_t::Stencil,
                                 Boundary_conditions_t>;

  class LB_boundary_handling {
  

  public:
    LB_boundary_handling(const walberla::BlockDataID &flag_field_id,
                         const walberla::BlockDataID &pdf_field_id,
                         const walberla::real_t top_velocity)
        : m_flag_field_id(flag_field_id), m_pdf_field_id(pdf_field_id),
          m_top_velocity(top_velocity) {}

    Boundary_handling_t *operator()(walberla::IBlock *const block) {

      Flag_field_t *flag_field = block->getData<Flag_field_t>(m_flag_field_id);
      Pdf_field_t *pdf_field = block->getData<Pdf_field_t>(m_pdf_field_id);

      const uint8_t fluid = flag_field->registerFlag(Fluid_flag);

      return new Boundary_handling_t(
          "boundary handling", flag_field, fluid,
          boost::tuples::make_tuple(
              No_slip_t("no slip", No_slip_flag, pdf_field),
              UBB_t("velocity bounce back", UBB_flag, pdf_field, m_top_velocity,
                    walberla::real_c(0), walberla::real_c(0))));
    }

  private:
    const walberla::BlockDataID m_flag_field_id;
    const walberla::BlockDataID m_pdf_field_id;

    const walberla::real_t m_top_velocity;
  };

public:
  LbWalberla(double viscosity, double agrid,
             const Vector3d &box_dimensions, const Vector3i &node_grid,
             double skin);

  void integrate();
  Vector3d get_node_velocity(const Vector3i node) const;
  void set_node_velocity(const Vector3i node, const Vector3d v);

  void add_force_at_pos(const Vector3d &position, const Vector3d &force);
  Vector3d get_velocity_at_pos(const Vector3d &position) const;
  double get_density_at_pos(const Vector3d &position);
  //  Vector3d get_stress_at_pos(const Vector3d& position);

  void set_node_as_boundary(const Vector3i &node);
  //  void set_node_as_boundary(const Vector3i &node,
  //                            const Vector3d &slip_velocity);
  void remove_node_from_boundary(const Vector3i &node);
  int get_node_is_boundary(const Vector3i &node);

  void print_vtk_velocity(char *filename);
  void print_vtk_density(char *filename);
  void print_vtk_boundary(char *filename);

  void set_external_force(const Vector3d &ext_force) {
    m_ext_force = ext_force;
  }
  Vector3d get_external_force() { return m_ext_force; }

  void set_viscosity(double viscosity);
  double get_viscosity();


  Vector3i get_grid_dimensions() {
    return Vector3i{{int(m_blocks->getXSize()), int(m_blocks->getYSize()), int(m_blocks->getZSize())}};
  };

  double get_grid_spacing() {
    return m_agrid;
  }

private:

  walberla::BlockDataID m_pdf_field_id;
  walberla::BlockDataID m_flag_field_id;
  walberla::BlockDataID m_force_field_id;
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

  std::shared_ptr<walberla::vtk::VTKOutput> create_fluid_field_vtk_writer(
      std::shared_ptr<walberla::blockforest::StructuredBlockForest> &blocks,
      const walberla::BlockDataID &pdf_field_id,
      const walberla::BlockDataID &flag_field_id);
};

void walberla_mpi_init();


#endif // LB_WALBERLA

#endif // LB_WALBERLA_H
