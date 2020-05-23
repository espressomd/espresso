#include "LBBoundary.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
namespace LBBoundaries {
Utils::Vector3d LBBoundary::get_force() const {
#ifdef LB_BOUNDARIES
  if (lattice_switch == ActiveLB::WALBERLA) {
#ifdef LB_WALBERLA
    auto const grid = lb_walberla()->get_grid_dimensions();
    auto const agrid = lb_lbfluid_get_agrid();
    Utils::Vector3d force{0, 0, 0};
    for (auto index_and_pos : lb_walberla()->node_indices_positions(true)) {
      // Convert to MD units
      auto const index = index_and_pos.first;
      auto const pos = index_and_pos.second * agrid;
      if (not shape().is_inside(pos)) {
             auto node_force_density =
                  lb_walberla()->get_node_boundary_force(index);
              if (node_force_density) {
                force += (*node_force_density);
              }

            }
    }         // loop over lb cells
    return boost::mpi::all_reduce(comm_cart, force,
                                  std::plus<Utils::Vector3d>());
#endif
  }
#endif
  throw std::runtime_error("LB Boundary code called with inactive LB");
}
} // namespace LBBoundaries
