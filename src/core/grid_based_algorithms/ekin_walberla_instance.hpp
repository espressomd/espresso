#ifndef ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP
#define ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP

#include "EKinWalberlaBase.hpp"
#include "communication.hpp"
#include "ekin_walberla_init.hpp"
#include "ekin_walberla_interface.hpp"

struct EKWalberlaParams {
  explicit EKWalberlaParams(double tau) : m_tau(tau) {}
  [[nodiscard]] double get_tau() const { return m_tau; };

private:
  double m_tau;
};

/** @brief Access the per-MPI-node EKin instance */
EKinWalberlaBase<double> *ekin_walberla();

/** @brief Access the EK Walberla parameters */
EKWalberlaParams *ek_walberla_params();

void mpi_init_ekin_walberla(double diffusion, double kT, double density,
                            double tau);
void mpi_destruct_ekin_walberla();

#endif // ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP
