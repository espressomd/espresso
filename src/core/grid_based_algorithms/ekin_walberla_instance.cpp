#include "config.hpp"

#ifdef EK_WALBERLA
#include "ekin_walberla_instance.hpp"

#include "walberla_blockforest.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"

#include "ekin_walberla_interface.hpp"

#include "EKinWalberlaBase.hpp"
#include "ekin_walberla_init.hpp"

#include <memory>

namespace {
EKinWalberlaBase<double> *ekin_walberla_instance = nullptr;
std::unique_ptr<EKWalberlaParams> ek_walberla_params_instance{nullptr};
} // namespace

EKinWalberlaBase<double> *ekin_walberla() {
  if (!ekin_walberla_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized EKinWalberla instance.");
  }
  return ekin_walberla_instance;
}

EKWalberlaParams *ek_walberla_params() {
  if (!ek_walberla_params_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized EKinWalberlaParams instance.");
  }
  return ek_walberla_params_instance.get();
}

void init_ekin_walberla_local(double diffusion, double kT, double density,
                              double tau) {
  // Exceptions need to be converted to runtime errors so they can be
  // handled from Python in a parallel simulation
  try {
    ekin_walberla_instance =
        new_ekin_walberla(get_walberla_blockforest(), diffusion, kT, density);
    ek_walberla_params_instance = std::make_unique<EKWalberlaParams>(tau);
  } catch (const std::exception &e) {
    runtimeErrorMsg() << "Error during Walberla initialization: " << e.what();
    ekin_walberla_instance = nullptr;
    ek_walberla_params_instance.reset();
  }
}
REGISTER_CALLBACK(init_ekin_walberla_local)

void destruct_ekin_walberla_local() {
  delete ekin_walberla_instance;
  ekin_walberla_instance = nullptr;
}
REGISTER_CALLBACK(destruct_ekin_walberla_local)

void mpi_init_ekin_walberla(double diffusion, double kT, double density,
                            double tau) {
  mpi_call_all(init_ekin_walberla_local, diffusion, kT, density, tau);
  mpi_set_ek_lattice_switch(EK::ActiveEK::WALBERLA);
}

void mpi_destruct_ekin_walberla() {
  Communication::mpiCallbacks().call_all(destruct_ekin_walberla_local);
  mpi_set_ek_lattice_switch(EK::ActiveEK::NONE);
}
#endif // EK_WALBERLA