#include "config.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "DipolarDirectSum_cuda.hpp"
#include "DipolarDirectSum.hpp"
#include "../forces.hpp"
#include "EspressoSystemInterface.hpp"
#include "forces.hpp"
#include "energy.hpp"




#ifdef DIPOLAR_DIRECT_SUM

void activate_dipolar_direct_sum_gpu()
{
if (dipolarDirectSum)
  free(dipolarDirectSum);

// also necessary on 1 CPU or GPU, does more than just broadcasting
mpi_bcast_coulomb_params();
dipolarDirectSum =new DipolarDirectSum(espressoSystemInterface);
forceActors.push_back(dipolarDirectSum);
energyActors.push_back(dipolarDirectSum);

coulomb.Dmethod = DIPOLAR_DS_GPU;

}

void deactivate_dipolar_direct_sum_gpu()
{
if (dipolarDirectSum)
{
  forceActors.remove(dipolarDirectSum);
  energyActors.remove(dipolarDirectSum);
  delete(dipolarDirectSum);
  dipolarDirectSum=nullptr;

}
}


DipolarDirectSum *dipolarDirectSum=0;

#endif

