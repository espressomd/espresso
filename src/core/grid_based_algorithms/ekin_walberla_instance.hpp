#ifndef ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP
#define ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP

#include "EKinWalberlaBase.hpp"
#include "communication.hpp"
#include "ekin_walberla_init.hpp"
#include "ekin_walberla_interface.hpp"

/** @brief Access the per-MPI-node EKin instance */
EKinWalberlaBase<double> *ekin_walberla();

void mpi_init_ekin_walberla(double diffusion, double kT, double density);
void mpi_destruct_ekin_walberla();

#endif // ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP
