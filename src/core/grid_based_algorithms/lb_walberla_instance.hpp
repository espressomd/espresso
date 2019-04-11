#ifndef GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP
#define GRID_BASED_ALGORITHMS_LBWALBERLA_INSTANCE_HPP

#include "LbWalberla.hpp"
#include "utils/Vector.hpp"

/** @brief Initialize Walberla's MPI manager */
void walberla_mpi_init();

/** @brief Access the per-MPI-node LbWalberla isntance */
const LbWalberla *lb_walberla();

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

#endif
