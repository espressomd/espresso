#include "DipolarDirectSum.hpp"

#include "EspressoSystemInterface.hpp"
#include "forces.hpp"
#include "energy.hpp"

#include "utils/make_unique.hpp"

#ifdef DIPOLAR_DIRECT_SUM

std::unique_ptr<DipolarDirectSum> dipolarDirectSum;

void activate_dipolar_direct_sum_gpu() {
  // also necessary on 1 CPU or GPU, does more than just broadcasting
  coulomb.Dmethod = DIPOLAR_DS_GPU;
  mpi_bcast_coulomb_params();

  dipolarDirectSum = Utils::make_unique<DipolarDirectSum>(espressoSystemInterface);
  forceActors.push_back(dipolarDirectSum.get());
  energyActors.push_back(dipolarDirectSum.get());
}

void deactivate_dipolar_direct_sum_gpu() {
  if (dipolarDirectSum) {
    forceActors.remove(dipolarDirectSum.get());
    energyActors.remove(dipolarDirectSum.get());
    dipolarDirectSum = nullptr;
  }
}

#endif
