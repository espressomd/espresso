#include "config.hpp"

#ifdef EK_WALBERLA
#include "ekin_walberla_instance.hpp"

#include "walberla_blockforest.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"

#include "ek_interface.hpp"

#include "EKinWalberlaBase.hpp"
#include "ekin_walberla_init.hpp"

#include <cassert>
#include <memory>

namespace {
std::vector<EKWalberlaInstance> ek_walberla_instances;
EKinWalberlaBase<double> *ekin_walberla_instance = nullptr;
std::unique_ptr<EKWalberlaParams> ek_walberla_params_instance{nullptr};
} // namespace

EKinWalberlaBase<double> *ekin_walberla(uint id) {
  return get_ek_instance_walberla(id).get_ek();
}

std::vector<EKWalberlaInstance> &get_eks_walberla() {
  return ek_walberla_instances;
}

const EKWalberlaInstance &get_ek_instance_walberla(uint id) {
  assert(!get_eks_walberla().empty());

  return get_eks_walberla().at(id);
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
    get_eks_walberla().emplace_back(
        new_ekin_walberla(get_walberla_blockforest(), diffusion, kT, density),
        tau);
  } catch (const std::exception &e) {
    runtimeErrorMsg() << "Error during Walberla initialization: " << e.what();
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