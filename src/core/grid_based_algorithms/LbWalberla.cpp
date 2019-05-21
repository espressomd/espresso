#include "config.hpp"
#ifdef LB_WALBERLA
#include "LbWalberla.hpp"
#include "communication.hpp"
#include "lb_walberla_instance.hpp"
#include "utils/Vector.hpp"
#include "utils/mpi/gatherv.hpp"

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
#include "field/distributors/DistributorCreators.h"
#include "field/interpolators/FieldInterpolatorCreators.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/UBB.h"
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

#include "boost/optional.hpp"

using namespace walberla;

inline Utils::Vector3d to_vector3d(const Vector3<real_t> v) {
  return Utils::Vector3d{v[0], v[1], v[2]};
}
Vector3<real_t> to_vector3(const Utils::Vector3d v) {
  return Vector3<real_t>{v[0], v[1], v[2]};
}

LbWalberla::LbWalberla(double viscosity, double density, double agrid,
                       double tau, const Utils::Vector3d &box_dimensions,
                       const Utils::Vector3i &node_grid, double skin) {

  m_skin = skin;
  m_agrid = agrid;
  m_tau = tau;
  m_density = density;
  m_velocity = {0, 0, 0};

  Utils::Vector3i grid_dimensions;
  for (int i = 0; i < 3; i++) {
    if (fabs(floor(box_dimensions[i] / agrid) * agrid - box_dimensions[i]) >
        std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error(
          "Box length not commensurate with agrid in direction " +
          std::to_string(i));
    }
    grid_dimensions[i] = int(std::round(box_dimensions[i] / agrid));
    if (grid_dimensions[i] % node_grid[i] != 0) {
      throw std::runtime_error(
          "LB grid dimensions and mpi node grid are not compatible.");
    }
  }
  m_grid_dimensions = grid_dimensions;

  m_blocks = blockforest::createUniformBlockGrid(
      uint_c(node_grid[0]), // blocks in x direction
      uint_c(node_grid[1]), // blocks in y direction
      uint_c(node_grid[2]), // blocks in z direction
      uint_c(grid_dimensions[0] /
             node_grid[0]), // number of cells per block in x direction
      uint_c(grid_dimensions[1] /
             node_grid[1]), // number of cells per block in y direction
      uint_c(grid_dimensions[2] /
             node_grid[2]), // number of cells per block in z direction
      real_c(1.0),          // Lattice constant
      uint_c(node_grid[0]), uint_c(node_grid[1]),
      uint_c(node_grid[2]), // cpus per direction
      true, true, true, true);

  m_force_field_id = field::addToStorage<vector_field_t>(
      m_blocks, "force field", math::Vector3<real_t>{0, 0, 0}, field::zyxf,
      uint_c(1));
  m_lattice_model = std::make_shared<Lattice_model_t>(Lattice_model_t(
      lbm::collision_model::TRT::constructWithMagicNumber(
          lbm::collision_model::omegaFromViscosity((real_t)viscosity)),
      force_model_t(m_force_field_id)));

  m_pdf_field_id =
      lbm::addPdfFieldToStorage(m_blocks, "pdf field", *m_lattice_model,
                                to_vector3(m_velocity), (real_t)m_density);

  m_flag_field_id =
      field::addFlagFieldToStorage<Flag_field_t>(m_blocks, "flag field");

  m_boundary_handling_id = m_blocks->addBlockData<Boundary_handling_t>(
      LB_boundary_handling(m_flag_field_id, m_pdf_field_id),
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

boost::optional<LbWalberla::BlockAndCell>
LbWalberla::get_block_and_cell(const Utils::Vector3i &node) const {
  // Get block and local cell
  Cell global_cell{uint_c(node[0]), uint_c(node[1]), uint_c(node[2])};
  auto block = m_blocks->getBlock(global_cell, 0);
  // Return if we don't have the cell
  if (!block)
    return {boost::none};

  // Transform coords to block local
  Cell local_cell;
  m_blocks->transformGlobalToBlockLocalCell(local_cell, *block, global_cell);
  return {{block, local_cell}};
}

bool LbWalberla::set_node_velocity_at_boundary(const Utils::Vector3i node,
                                               const Utils::Vector3d v) {
  auto bc = get_block_and_cell(node);
  // Return if we don't have the cell.
  if (!bc)
    return false;

  const UBB_t::Velocity velocity(real_c(v[0]), real_c(v[1]), real_c(v[2]));

  Boundary_handling_t *boundary_handling =
      (*bc).block->getData<Boundary_handling_t>(m_boundary_handling_id);
  boundary_handling->forceBoundary(UBB_flag, bc->cell[0], bc->cell[1],
                                   bc->cell[2], velocity);
  return true;
}

boost::optional<Utils::Vector3d>
LbWalberla::get_node_velocity_at_boundary(const Utils::Vector3i &node) const {
  boost::optional<Utils::Vector3d> res;
  auto bc = get_block_and_cell(node);
  // return if we don't have the cell
  if (!bc)
    return res;
  const Boundary_handling_t *boundary_handling =
      (*bc).block->getData<Boundary_handling_t>(m_boundary_handling_id);
  walberla::boundary::BoundaryUID uid =
      boundary_handling->getBoundaryUID(UBB_flag);

  if (!boundary_handling->isBoundary((*bc).cell))
    return res;

  Utils::Vector3d v =
      to_vector3d(boundary_handling->getBoundaryCondition<UBB_t>(uid).getValue(
          (*bc).cell[0], (*bc).cell[1], (*bc).cell[2]));
  res = {v};
  return res;
}

bool LbWalberla::remove_node_from_boundary(const Utils::Vector3i &node) {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return false;
  Boundary_handling_t *boundary_handling =
      (*bc).block->getData<Boundary_handling_t>(m_boundary_handling_id);
  boundary_handling->removeBoundary((*bc).cell[0], (*bc).cell[1],
                                    (*bc).cell[2]);
  return true;
}

boost::optional<bool>
LbWalberla::get_node_is_boundary(const Utils::Vector3i &node) const {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return {boost::none};

  Boundary_handling_t *boundary_handling =
      (*bc).block->getData<Boundary_handling_t>(m_boundary_handling_id);
  return {boundary_handling->isBoundary((*bc).cell)};
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

boost::optional<Utils::Vector3d>
LbWalberla::get_node_velocity(const Utils::Vector3i node) const {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return {boost::none};
  // Get pdf field
  //     auto const& pdf_field = block->getData<Pdf_field_t>(m_pdf_field_id);
  auto const &vel_adaptor =
      (*bc).block->getData<VelocityAdaptor>(m_velocity_adaptor_id);
  return {to_vector3d(vel_adaptor->get((*bc).cell))};
}

boost::optional<Utils::Vector3d>
LbWalberla::get_velocity_at_pos(const Utils::Vector3d &pos) const {
  auto block = m_blocks->getBlock(to_vector3(pos));
  if (!block)
    return {boost::none};

  auto *velocity_interpolator = block->getData<VectorFieldAdaptorInterpolator>(
      m_velocity_interpolator_id);
  Vector3<real_t> v;
  velocity_interpolator->get(to_vector3(pos), &v);
  return {to_vector3d(v)};
}

bool LbWalberla::set_node_velocity(const Utils::Vector3i &node,
                                   const Utils::Vector3d v) {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return false;

  auto pdf_field = (*bc).block->getData<Pdf_field_t>(m_pdf_field_id);
  const real_t density = pdf_field->getDensity((*bc).cell);

  pdf_field->setDensityAndVelocity((*bc).cell,
                                   Vector3<double>{v[0], v[1], v[2]}, density);
  return true;
}

bool LbWalberla::set_node_density(const Utils::Vector3i node,
                                  const double density) {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return false;

  auto pdf_field = (*bc).block->getData<Pdf_field_t>(m_pdf_field_id);
  auto const &vel_adaptor =
      (*bc).block->getData<VelocityAdaptor>(m_velocity_adaptor_id);
  Vector3<double> v = vel_adaptor->get((*bc).cell);

  pdf_field->setDensityAndVelocity((*bc).cell,
                                   Vector3<double>{v[0], v[1], v[2]}, density);

  return true;
}

boost::optional<double>
LbWalberla::get_node_density(const Utils::Vector3i node) const {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return {boost::none};

  auto pdf_field = (*bc).block->getData<Pdf_field_t>(m_pdf_field_id);

  return {pdf_field->getDensity((*bc).cell)};
}

bool LbWalberla::set_node_pop(const Utils::Vector3i &node,
                              const Utils::Vector19d pop) {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return false;

  auto pdf_field = (*bc).block->getData<Pdf_field_t>(m_pdf_field_id);
  auto &xyz0 = (*pdf_field)(node[0], node[1], node[2], 0);

  for (int i = 0; i < 19; i++) {
    pdf_field->getF(&xyz0, es_pop_index_to_walberla_pop_index[i]) = pop[i];
  }

  return true;
}

boost::optional<Utils::Vector19d>
LbWalberla::get_node_pop(const Utils::Vector3i node) const {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return {boost::none};

  auto pdf_field = (*bc).block->getData<Pdf_field_t>(m_pdf_field_id);
  const auto &xyz0 = (*pdf_field)(node[0], node[1], node[2], 0);

  Utils::Vector19d pop;

  for (int i = 0; i < 19; i++) {
    pop[i] = pdf_field->getF(&xyz0, es_pop_index_to_walberla_pop_index[i]);
  }

  return {pop};
}

void LbWalberla::set_viscosity(double viscosity) {
  m_lattice_model->collisionModel().resetWithMagicNumber(
      lbm::collision_model::omegaFromViscosity((real_t)viscosity));
}

double LbWalberla::get_viscosity() {
  return m_lattice_model->collisionModel().viscosity();
}
bool LbWalberla::node_in_local_domain(const Utils::Vector3i &node) const {
  auto block = m_blocks->getBlock(
      Cell{uint_c(node[0]), uint_c(node[1]), uint_c(node[2])}, 0);
  return (block != nullptr);
}

bool LbWalberla::pos_in_local_domain(const Utils::Vector3d &pos) const {
  auto block =
      m_blocks->getBlock(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]));
  return (block != nullptr);
}

#endif
