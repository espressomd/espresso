#ifndef ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP
#define ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP

#include "config.hpp"

#ifdef EK_WALBERLA
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

struct EKWalberlaInstance {
  EKWalberlaInstance(EKinWalberlaBase<double> *ekptr, double tau)
      : m_tau(tau), m_ekptr(ekptr) {}
  [[nodiscard]] double get_tau() const { return m_tau; };
  [[nodiscard]] EKinWalberlaBase<double> *get_ek() const { return m_ekptr; };

private:
  double m_tau;
  EKinWalberlaBase<double> *m_ekptr;
};

/** @brief Access the per-MPI-node EKin instance */
EKinWalberlaBase<double> *ekin_walberla(uint id);

std::vector<EKWalberlaInstance> &get_eks_walberla();

const EKWalberlaInstance &get_ek_instance_walberla(uint id);

/** @brief Access the EK Walberla parameters */
EKWalberlaParams *ek_walberla_params();

void mpi_init_ekin_walberla(double diffusion, double kT, double density,
                            double tau);
void mpi_destruct_ekin_walberla();
#endif // EK_WALBERLA
#endif // ESPRESSO_EKIN_WALBERLA_INSTANCE_HPP
