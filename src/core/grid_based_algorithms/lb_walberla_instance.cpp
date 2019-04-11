#ifderf LB_WALBERLA
#include "LbWalberla.hpp"
#include "communication.hpp"
#include "utils/Vector.hpp"
#include <memory>

void walberla_mpi_init() {
  int argc = 0;
  char **argv = NULL;
  static mpi::Environment m_env = mpi::Environment(argc, argv);
}
namespace {
std::unique_ptr<LbWalberla> lb_walberla_instance = nullptr;
}

const LbWalberla *lb_walberla() {
  if (!lb_walberla_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized LbWalberla instance.");
  }
  return lb_walberla_instance.get();
}

void init_lb_walberla(double viscosity, double agrid,
                      const Vector3d &box_dimensions, const Vector3i &node_grid,
                      double skin) {
  lb_walberla_instance = std::make_unique<LbWalberla>(
      LbWalberla{viscosity, agrid, box_dimensions, node_grid, skin});
}
REGISTER_CALLBACK(init_lb_walberla);

void destruct_lb_walberla() { lb_walberla_instance.reset(nullptr); }
REGISTER_CALLBACK(destruct_lb_walberla);

#endif
