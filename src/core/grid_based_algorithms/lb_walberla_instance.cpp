#include "config.hpp"
#include "errorhandling.hpp"

#ifdef LB_WALBERLA
#include "LbWalberla.hpp"

#include "communication.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "lb_interface.hpp"

#include "utils/Vector.hpp"

#include <memory>

#include "core/mpi/Environment.h"

#ifdef LB_WALBERLA
void walberla_mpi_init() {
  int argc = 0;
  char **argv = NULL;
  static walberla::mpi::Environment m_env =
      walberla::mpi::Environment(argc, argv);
}

namespace {
std::unique_ptr<LbWalberla> lb_walberla_instance = nullptr;
}

LbWalberla *lb_walberla() {
  if (!lb_walberla_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized LbWalberla instance.");
  }
  return lb_walberla_instance.get();
}

void init_lb_walberla(double viscosity, double density, double agrid,
                      double tau, const Utils::Vector3d &box_dimensions,
                      const Utils::Vector3i &node_grid, double skin) {
  // Exceptions need to be converted to runtime erros so they can be
  // handled from Python in a parallel simulation
  try {
    lb_walberla_instance = std::make_unique<LbWalberla>(LbWalberla{
        viscosity, density, agrid, tau, box_dimensions, node_grid, skin});
  } catch (const std::exception &e) {
    runtimeErrorMsg() << "Error during Walberla initialization: " << e.what();
    lb_walberla_instance.reset(nullptr);
  }
}
REGISTER_CALLBACK(init_lb_walberla)

void destruct_lb_walberla() { lb_walberla_instance.reset(nullptr); }
REGISTER_CALLBACK(destruct_lb_walberla)

void mpi_init_lb_walberla(double viscosity, double density, double agrid,
                          double tau) {
  Communication::mpiCallbacks().call_all(init_lb_walberla, viscosity,
                                         density * pow(agrid, 3), agrid, tau,
                                         box_geo.length(), node_grid, skin);
  if (lb_walberla_instance) {
    lb_lbfluid_set_lattice_switch(ActiveLB::WALBERLA);
     lb_lbfluid_sanity_checks();
  }
}

void mpi_destruct_lb_walberla() {
  lb_lbfluid_set_lattice_switch(ActiveLB::NONE);
  Communication::mpiCallbacks().call_all(destruct_lb_walberla);
}
#endif
#endif
