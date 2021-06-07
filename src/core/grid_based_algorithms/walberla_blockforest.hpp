#ifndef ESPRESSO_WALBERLA_BLOCKFOREST_HPP
#define ESPRESSO_WALBERLA_BLOCKFOREST_HPP

#include "WalberlaBlockForest.hpp"

/** @brief Access the per-MPI-node Walberla instance */
std::shared_ptr<walberla::WalberlaBlockForest> walberla_blockforest;

void mpi_init_walberla_blockforest(const Utils::Vector3d &box_size,
                                   double agrid, int n_ghost_layers);

#endif // ESPRESSO_WALBERLA_BLOCKFOREST_HPP
