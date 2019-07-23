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
#include "core/mpi/MPITextFile.h"
#include "core/mpi/Reduce.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "field/distributors/DistributorCreators.h"
#include "field/interpolators/FieldInterpolatorCreators.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/UBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <memory>

#include "boost/optional.hpp"

using namespace walberla;

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
      printf("Grid dimension: %d, node grid %d\n",grid_dimensions[i],node_grid[i]);
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
      n_ghost_layers());
  m_force_field_from_md_id = field::addToStorage<vector_field_t>(
      m_blocks, "force field", math::Vector3<real_t>{0, 0, 0}, field::zyxf,
      n_ghost_layers());
  m_lattice_model = std::make_shared<Lattice_model_t>(Lattice_model_t(
      lbm::collision_model::TRT::constructWithMagicNumber(
          lbm::collision_model::omegaFromViscosity((real_t)viscosity)),
      force_model_t(m_force_field_id)));

  m_pdf_field_id = lbm::addPdfFieldToStorage(
      m_blocks, "pdf field", *m_lattice_model, to_vector3(m_velocity),
      (real_t)m_density, n_ghost_layers());

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
      std::make_shared<field::communication::PackInfo<Pdf_field_t>>(
          m_pdf_field_id));
  m_reset_force = std::make_shared<
      ResetForce<Pdf_field_t, vector_field_t, Boundary_handling_t>>(
      m_pdf_field_id, m_force_field_id, m_force_field_from_md_id,
      m_boundary_handling_id);
  m_time_loop->add() << // timeloop::BeforeFunction(communication, "communication")
                     timeloop::Sweep(Boundary_handling_t::getBlockSweep(
                                            m_boundary_handling_id),
                                        "boundary handling");
  m_time_loop->add() << timeloop::Sweep(makeSharedSweep(m_reset_force),
                                        "Reset force fields");
  m_time_loop->add() << timeloop::AfterFunction(communication, "communication")
     << timeloop::Sweep(
      domain_decomposition::makeSharedSweep(
          lbm::makeCellwiseSweep<Lattice_model_t, Flag_field_t>(
              m_pdf_field_id, m_flag_field_id, Fluid_flag)),
      "LB stream & collide");

  m_force_distributor_id =
      field::addDistributor<Vector_field_distributor_t, Flag_field_t>(
          m_blocks, m_force_field_from_md_id, m_flag_field_id, Fluid_flag);

  m_velocity_adaptor_id = field::addFieldAdaptor<VelocityAdaptor>(
      m_blocks, m_pdf_field_id, "velocity adaptor");

  m_velocity_interpolator_id =
      field::addFieldInterpolator<VectorFieldAdaptorInterpolator, Flag_field_t>(
          m_blocks, m_velocity_adaptor_id, m_flag_field_id, Fluid_flag);

  communication();
}

Utils::Vector3d LbWalberla::get_momentum() const {
  Vector3<real_t> mom;
  for (auto block_it = m_blocks->begin(); block_it != m_blocks->end();
       ++block_it) {
    auto pdf_field = block_it->getData<Pdf_field_t>(m_pdf_field_id);
    Vector3<real_t> local_v;
    WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
      double local_dens = pdf_field->getDensityAndVelocity(local_v, x, y, z);
      mom += local_dens * local_v;
    });
  }
  return to_vector3d(mom);
}

void LbWalberla::integrate() { m_time_loop->singleStep(); }

boost::optional<LbWalberla::BlockAndCell>
LbWalberla::get_block_and_cell(const Utils::Vector3i &node, bool consider_ghost_layers) const {
  // Get block and local cell
  Cell global_cell{uint_c(node[0]), uint_c(node[1]), uint_c(node[2])};
  auto block = m_blocks->getBlock(global_cell, 0);
  // Return if we don't have the cell
  if (consider_ghost_layers and !block) {
    // Try to find a block which has the cell as ghost layer
    for (auto b = m_blocks->begin(); b!=m_blocks->end(); ++b) {
      if (b->getAABB().getExtended(m_skin/m_agrid).contains(real_c(node[0]),real_c(node[1]),real_c(node[2]))) {
        block=&(*b);
        break;
      }
    }
  }
  if (!block)
    return {boost::none};

  // Transform coords to block local
  Cell local_cell;
  m_blocks->transformGlobalToBlockLocalCell(local_cell, *block, global_cell);
  return {{block, local_cell}};
}

IBlock* 
LbWalberla::get_block(const Utils::Vector3d &pos, bool consider_ghost_layers) const {
  // Get block and local cell
  auto block = m_blocks->getBlock(real_c(pos[0]), real_c(pos[1]), real_c(pos[2]));
  if (consider_ghost_layers and !block) {
    // Try to find a block which has the cell as ghost layer
    for (auto b = m_blocks->begin(); b!=m_blocks->end(); ++b) {
      auto ab =b->getAABB().getExtended(1.0 *n_ghost_layers());
      printf("%d: %g %g %g, %g %g %g\n",this_node, ab.min()[0],ab.min()[1],ab.min()[2],ab.max()[0],ab.max()[1],ab.max()[2]);
      printf("%d: pos %g %g %g\n",this_node, pos[0],pos[1],pos[2]);
      if (ab.contains(pos[0],pos[1],pos[2])) {
        block=&(*b);
        printf("%d: block selected\n",this_node);
        break;
      }
    }
  }

  return block;
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

bool LbWalberla::add_force_at_pos(const Utils::Vector3d &pos,
                                  const Utils::Vector3d &force) {
  auto block = get_block(pos,true);
  if (!block)
    return false;
  auto *force_distributor =
      block->getData<Vector_field_distributor_t>(m_force_distributor_id);
  auto f = to_vector3(force);
  force_distributor->distribute(to_vector3(pos), &f);
  return true;
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
  auto block = get_block(pos,true);
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

boost::optional<Utils::Vector6d>
LbWalberla::get_node_pressure_tensor(const Utils::Vector3i node) const {
  auto bc = get_block_and_cell(node);
  if (!bc)
    return {boost::none};
  auto pdf_field = (*bc).block->getData<Pdf_field_t>(m_pdf_field_id);

  return to_vector6d(pdf_field->getPressureTensor((*bc).cell));
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
