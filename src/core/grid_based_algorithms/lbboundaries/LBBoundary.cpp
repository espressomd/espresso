#include "LBBoundary.hpp"
#include "config.hpp"
#ifdef LB_BOUNDARIES
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_walberla_instance.hpp"
namespace LBBoundaries {
Utils::Vector3d LBBoundary::get_force() const {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
  if (lattice_switch == ActiveLB::WALBERLA) {
    auto const grid = lb_walberla()->get_grid_dimensions();
    auto const agrid = lb_lbfluid_get_agrid();
    Utils::Vector3d force_density{0, 0, 0};
    for (auto index_and_pos : lb_walberla()->global_node_indices_positions()) {
      // Convert to MD units
      auto const index = index_and_pos.first;
      auto const pos = index_and_pos.second * agrid;
      if (calc_dist(pos) <= 0.)
        for (int dx : {-1, 0, 1})
          for (int dy : {-1, 0, 1})
            for (int dz : {-1, 0, 1}) {
              Utils::Vector3i shifted_index =
                  index +
                  Utils::Vector3i{dx * grid[0], dy * grid[1], dz * grid[2]};
              auto node_force_density =
                  lb_walberla()->get_node_boundary_force(shifted_index);
              if (node_force_density) {
                force_density += (*node_force_density);
              }
            } // Loop over cells
    }         // loop over lb cells
    return boost::mpi::all_reduce(comm_cart, force_density,
                                  std::plus<Utils::Vector3d>()) *
           lb_lbfluid_get_density();
  } else {
    if (this_node == 0) {
      return lbboundary_get_force(this);
    } else
      return {};
  }
#else
  throw std::runtime_error("Needs LB_BOUNDARIES or LB_BOUNDARIES_GPU.");
#endif
}
} // namespace LBBoundaries
#endif
