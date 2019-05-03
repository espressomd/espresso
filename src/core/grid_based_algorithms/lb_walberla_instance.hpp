#ifndef GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP
#define GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP
#include "config.hpp"

#ifdef LB_WALBERLA
#include "LbWalberla.hpp"
#include "communication.hpp"
#include "utils/Vector.hpp"
#include "utils/mpi/gatherv.hpp"

/** @brief Initialize Walberla's MPI manager */
void walberla_mpi_init();

/** @brief Access the per-MPI-node LbWalberla isntance */
LbWalberla *lb_walberla();

/** @brief Create the LbWalberal instance and sets the lattice dwitch to WALBERLA
 *
 *  @param viscosity Fluid viscosity
 *  @param density Fluiddensity
 *  @param agrid  Size of one lb cell
 *  @param tau    LB time step
 */
void mpi_init_lb_walberla(double viscosity, double density, double agrid, double tau);

/** @brief Destruct the LbWalberal instance and set lattice switch to NONE */
void mpi_destruct_lb_walberla();

/** @brief Get local fluid properties of the Walberla instance on the head node
 * @param getter  reference to the getter function, e.g.,
 * &LbWalberla::get_node_density The getter function has to take a single
 * argument *e.g., node index) and return a boost::optional indicating whether
 * the result is locally available.
 * @param arg     Single argument which is passed to the getter function
 * This function has to be called on all MPI nodes. It calls the getter
 * on the active LbWalberla instance as given by lb_walberla().
 *
 * @ret On the head node, the return value of the node-local getter is returned.
 * On other nodes, the return value is undefined.
 */
template <typename Getter, typename T>
auto lb_walberla_get_on_head_node(Getter getter_function, T const &arg)
    -> decltype(auto) {
  auto result = (lb_walberla()->*getter_function)(arg);

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
  typedef typename decltype(result)::value_type ResultType;
  ResultType global_result{};
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

#endif
