#ifdef LB_WALBERLA
#include "LbWalberla.hpp"
#include "utils/Vector.hpp"


#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimpleUBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <memory>

using namespace walberla;
int argc = 0;
char **argv = NULL;
mpi::Environment m_env = mpi::Environment(argc, argv);

inline Vector3d to_vector3d(const Vector3<real_t> v) {
  return Vector3d{v[0], v[1], v[2]};
}
Vector3<real_t> to_vector3(const Vector3d v) {
  return Vector3<real_t>{v[0], v[1], v[2]};
}

LbWalberla::LbWalberla(double viscosity, double agrid,
                       const Vector3d &box_dimensions,
                       const Vector3i &node_grid, double skin) {
  
  m_skin = skin;
  
  const Vector3i grid_dimensions = box_imensions/agrid;

  // Probably needs some args
  int argc = 0;
  m_blocks = blockforest::createUniformBlockGrid(
      node_grid[0], // blocks in x direction
      node_grid[1], // blocks in y direction
      node_grid[2], // blocks in z direction
      grid_dimensions[0] /
          node_grid[0], // number of cells per block in x direction
      grid_dimensions[1] /
          node_grid[1], // number of cells per block in y direction
      grid_dimensions[2] /
          node_grid[2], // number of cells per block in z direction
      1,                // Lattice constant
      false);           // more than one block can be on the same process

  m_force_field_id = field::addToStorage<vector_field_t>(
      m_blocks, "force field", math::Vector3<real_t>{0, 0, 0}, field::zyxf,
      uint_c(1));
  m_lattice_model = std::make_shared<Lattice_model_t>(
      lbm::collision_model::SRT(
          lbm::collision_model::omegaFromViscosity((real_t)viscosity)),
      force_model_t(m_force_field_id));

  m_pdf_field_id =
      lbm::addPdfFieldToStorage(m_blocks, "pdf field", *m_lattice_model);

  m_flag_field_id =
      field::addFlagFieldToStorage<Flag_field_t>(m_blocks, "flag field");

  m_boundary_handling_id = m_blocks->addBlockData<Boundary_handling_t>(
      LB_boundary_handling(m_flag_field_id, m_pdf_field_id, 0.0),
      "boundary handling");

  empty_flags_for_boundary(m_blocks, m_boundary_handling_id);

  // 1 timestep each time the integrate function is called
  m_time_loop =
      std::make_shared<timeloop::SweepTimeloop>(m_blocks->getBlockStorage(), 1);

  // sets up the communication
  blockforest::communication::UniformBufferedScheme<
      Lattice_model_t::CommunicationStencil>
      communication(m_blocks);
  communication.addPackInfo(
      std::make_shared<lbm::PdfFieldPackInfo<Lattice_model_t>>(m_pdf_field_id));

  m_time_loop->add() << timeloop::BeforeFunction(communication, "communication")
                     << timeloop::Sweep(Boundary_handling_t::getBlockSweep(
                                            m_boundary_handling_id),
                                        "boundary handling");
  m_time_loop->add() << timeloop::Sweep(
      domain_decomposition::makeSharedSweep(
          lbm::makeCellwiseSweep<Lattice_model_t, Flag_field_t>(
              m_pdf_field_id, m_flag_field_id, Fluid_flag)),
      "LB stream & collide");

  m_force_distributor_id =
      field::addDistributor<Vector_field_distributor_t, Flag_field_t>(
          m_blocks, m_force_field_id, m_flag_field_id, Fluid_flag);

  m_velocity_adaptor_id = field::addFieldAdaptor<VelocityAdaptor>(
      m_blocks, m_pdf_field_id, "velocity adaptor");

  m_velocity_interpolator_id =
      field::addFieldInterpolator<VectorFieldAdaptorInterpolator, Flag_field_t>(
          m_blocks, m_velocity_adaptor_id, m_flag_field_id, Fluid_flag);
}

void LbWalberla::print_vtk_density(char *filename) {
  // TODO: set the output file name to filename
  auto pdf_field_vtk_writer =
      create_fluid_field_vtk_writer(m_blocks, m_pdf_field_id, m_flag_field_id);
  vtk::VTKOutput::Write output = vtk::writeFiles(pdf_field_vtk_writer);
  output.forceWriteNextStep();
  output.execute();
}

void LbWalberla::integrate() { m_time_loop->run(); }

void LbWalberla::set_node_as_boundary(const Vector3i &node) {
  // TODO: Variable boundary flag and check length scale
  auto block = m_blocks->getBlock(node[0], node[1], node[2]);
  if (block != nullptr) {
    Boundary_handling_t *boundary_handling =
        block->getData<Boundary_handling_t>(m_boundary_handling_id);
    boundary_handling->forceBoundary(No_slip_flag, node[0], node[1], node[2]);
  }
}

void LbWalberla::remove_node_from_boundary(const Vector3i &node) {
  auto block = m_blocks->getBlock(node[0], node[1], node[2]);
  if (block != nullptr) {
    Boundary_handling_t *boundary_handling =
        block->getData<Boundary_handling_t>(m_boundary_handling_id);
    boundary_handling->removeBoundary(node[0], node[1], node[2]);
  }
}

int LbWalberla::get_node_is_boundary(const Vector3i &node) {
  auto block = m_blocks->getBlock(node[0], node[1], node[2]);
  if (block == nullptr)
    return -1;
  Boundary_handling_t *boundary_handling =
      block->getData<Boundary_handling_t>(m_boundary_handling_id);
  return boundary_handling->isBoundary(node[0], node[1], node[2]);
}

std::shared_ptr<walberla::vtk::VTKOutput>
LbWalberla::create_fluid_field_vtk_writer(
    std::shared_ptr<walberla::blockforest::StructuredBlockForest> &blocks,
    const walberla::BlockDataID &pdf_field_id,
    const walberla::BlockDataID &flag_field_id) {
  using namespace walberla;
  auto pdf_field_vtk_writer = vtk::createVTKOutput_BlockData(
      blocks, "fluid_field", uint_c(1), uint_c(1), true);

  blockforest::communication::UniformBufferedScheme<stencil::D3Q27>
      pdf_ghost_layer_sync(blocks);
  pdf_ghost_layer_sync.addPackInfo(
      make_shared<field::communication::PackInfo<Pdf_field_t>>(pdf_field_id));
  pdf_field_vtk_writer->addBeforeFunction(pdf_ghost_layer_sync);

  field::FlagFieldCellFilter<Flag_field_t> fluid_filter(flag_field_id);
  fluid_filter.addFlag(Fluid_flag);
  pdf_field_vtk_writer->addCellInclusionFilter(fluid_filter);

  auto velocity_writer =
      make_shared<lbm::VelocityVTKWriter<Lattice_model_t, float>>(
          pdf_field_id, "VelocityFromPDF");
  auto density_writer =
      make_shared<lbm::DensityVTKWriter<Lattice_model_t, float>>(
          pdf_field_id, "DensityFromPDF");
  pdf_field_vtk_writer->addCellDataWriter(velocity_writer);
  pdf_field_vtk_writer->addCellDataWriter(density_writer);

  return pdf_field_vtk_writer;
}

Vector3d LbWalberla::get_node_velocity(const Vector3i node) const {
  Cell global_cell{node[0], node[1], node[2]};
  // Get block which has the cell
  const IBlock *block = m_blocks->getBlock(global_cell, 0);

  // Transform coords to block local
  Cell local_cell;
  m_blocks->transformGlobalToBlockLocalCell(local_cell, *block, global_cell);

  // Get pdf field
  //     auto const& pdf_field = block->getData<Pdf_field_t>(m_pdf_field_id);
  auto const &vel_adaptor =
      block->getData<VelocityAdaptor>(m_velocity_adaptor_id);
  return to_vector3d(vel_adaptor->get(local_cell));
}

Vector3d LbWalberla::get_velocity_at_pos(const Vector3d &pos) const {
  auto block = m_blocks->getBlock(to_vector3(pos));
  auto *velocity_interpolator = block->getData<VectorFieldAdaptorInterpolator>(
      m_velocity_interpolator_id);
  Vector3<real_t> v;
  velocity_interpolator->get(to_vector3(pos), &v);
  return to_vector3d(v);
}

void LbWalberla::set_node_velocity(const Vector3i node, const Vector3d v) {
  Cell global_cell{node[0], node[1], node[2]};
  // Get block which has the cell
  auto block = m_blocks->getBlock(global_cell, 0);

  // Transform coords to block local
  Cell local_cell;
  m_blocks->transformGlobalToBlockLocalCell(local_cell, *block, global_cell);

  // Get pdf field
  auto pdf_field = block->getData<Pdf_field_t>(m_pdf_field_id);
  const real_t density = pdf_field->getDensity(local_cell);

  pdf_field->setDensityAndVelocity(local_cell,
                                   Vector3<double>{v[0], v[1], v[2]}, density);
}

void LbWalberla::set_viscosity(double viscosity) {
  m_lattice_model->collisionModel().reset(
      lbm::collision_model::omegaFromViscosity((real_t)viscosity));
}

double LbWalberla::get_viscosity() {
  return m_lattice_model->collisionModel().viscosity();
}

#endif LB_WALBERLA
