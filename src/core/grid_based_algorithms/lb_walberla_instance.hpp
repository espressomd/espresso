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

/** @brief Create the LbWalberal instance and sets the lattice dwitch to
 * WALBERLA
 *
 *  @param viscosity Fluid viscosity
 *  @param density Fluiddensity
 *  @param agrid  Size of one lb cell
 *  @param tau    LB time step
 */
void mpi_init_lb_walberla(double viscosity, double density, double agrid,
                          double tau);

/** @brief Destruct the LbWalberal instance and set lattice switch to NONE */
void mpi_destruct_lb_walberla();

#endif

#endif
