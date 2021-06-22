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
} // namespace

EKinWalberlaBase<double> *ekin_walberla() {
  if (!ekin_walberla_instance) {
    throw std::runtime_error(
        "Attempted access to uninitialized EKinWalberla instance.");
  }
  return ekin_walberla_instance;
}

void init_ekin_walberla_local(double diffusion, double kT, double density) {
  // Exceptions need to be converted to runtime errors so they can be
  // handled from Python in a parallel simulation
  try {
    ekin_walberla_instance =
        new_ekin_walberla(get_walberla_blockforest(), diffusion, kT, density);
  } catch (const std::exception &e) {
    runtimeErrorMsg() << "Error during Walberla initialization: " << e.what();
    ekin_walberla_instance = nullptr;
  }
}
REGISTER_CALLBACK(init_ekin_walberla_local)

void destruct_ekin_walberla_local() {
  delete ekin_walberla_instance;
  ekin_walberla_instance = nullptr;
}
REGISTER_CALLBACK(destruct_ekin_walberla_local)

void mpi_init_ekin_walberla(double diffusion, double kT, double density) {
  mpi_call_all(init_ekin_walberla_local, diffusion, kT, density);
  mpi_set_ek_lattice_switch(EK::ActiveEK::WALBERLA);
}

void mpi_destruct_ekin_walberla() {
  Communication::mpiCallbacks().call_all(destruct_ekin_walberla_local);
  mpi_set_ek_lattice_switch(EK::ActiveEK::NONE);
}