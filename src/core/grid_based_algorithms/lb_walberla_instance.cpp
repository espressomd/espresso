#include "config.hpp"

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
  lb_walberla_instance = std::make_unique<LbWalberla>(LbWalberla{
      viscosity, density, agrid, tau, box_dimensions, node_grid, skin});
}
REGISTER_CALLBACK(init_lb_walberla)

void destruct_lb_walberla() { lb_walberla_instance.reset(nullptr); }
REGISTER_CALLBACK(destruct_lb_walberla)

void mpi_init_lb_walberla(double viscosity, double density, double agrid,
                          double tau) {
  Communication::mpiCallbacks().call_all(init_lb_walberla, viscosity, density,
                                         agrid, tau, box_l, node_grid, skin);
  lb_lbfluid_set_lattice_switch(ActiveLB::WALBERLA);
}

void mpi_destruct_lb_walberla() {
  Communication::mpiCallbacks().call_all(destruct_lb_walberla);
  lb_lbfluid_set_lattice_switch(ActiveLB::NONE);
}
#endif
#endif
