#ifndef GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP
#define GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP

#include "LbWalberla.hpp"
#include "communication.hpp"
#include "utils/Vector.hpp"
#include "utils/mpi/gatherv.hpp"

/** @brief Initialize Walberla's MPI manager */
void walberla_mpi_init();

/** @brief Access the per-MPI-node LbWalberla isntance */
LbWalberla *lb_walberla();

/** @brief Create the per-MPI-rank LbWalberal instance
 *
 *  @param viscosity Fluid viscosity
 *  @param agrid  Size of one lb cell
 *  @param box_dimensions Dimensions of the global simulation box
 *  @param node_grid  Dimensions of the MPI node grid
 *  @param skin Distance beyond the node boundary a particle can reach
 */
void init_lb_walberla(double viscosity, double agrid,
                      const Vector3d &box_dimensions, const Vector3i &node_grid,
                      double skin);

/** @brief Destruct the per-MPI-rank LbWalberal instance */
void destruct_lb_walberla();

template <typename Getter, typename V>
auto collector_function(Getter getter_function, V position) -> decltype(auto) {
  auto result = (lb_walberla()->*getter_function)(position);

#ifdef ADDITIONAL_CHECKS
  int hits = 0;
  if (result) {
    hits = 1;
  }
  mpi_allreduce(&hits, &hits, 1, MPI_INT, MPI_SUM, comm_cart);
  if (hits > 1)
    throw std::runtime_error("More than one process found local data.");
  if (hits < 1)
    throw std::runtime_error("Local data cannot be found.");
#endif
  typedef typename decltype(result)::value_type T;
  T global_result{};
  if (comm_cart.rank() == 0) {
    if (!result) {
      comm_cart.recv(boost::mpi::any_source, 222, global_result);
    } else {
      global_result = result.get();
    }
    return global_result;
  } else {
    if (result) {
      global_result = result.get();
      comm_cart.send(0, 222, global_result);
    }
  }
  return global_result;
}

#endif
