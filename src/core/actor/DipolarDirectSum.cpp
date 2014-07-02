#include "config.hpp"
#include "grid.hpp"
#include "DipolarDirectSum_cuda.hpp"
#include "DipolarDirectSum.hpp"
#include "../forces.hpp"
#include "EspressoSystemInterface.hpp"
#include "forces.hpp"

#ifdef DIPOLAR_DIRECT_SUM




void activate_dipolar_direct_sum_gpu()
{
if (dipolarDirectSum)
  free(dipolarDirectSum);

dipolarDirectSum =new DipolarDirectSum(espressoSystemInterface);
potentials.push_back(dipolarDirectSum);

coulomb.Dmethod = DIPOLAR_DS_GPU;

}

void deactivate_dipolar_direct_sum_gpu()
{
if (dipolarDirectSum)
{
  potentials.remove(dipolarDirectSum);
  free(dipolarDirectSum);
  coulomb.Dmethod = DIPOLAR_NONE;
 }
}



DipolarDirectSum *dipolarDirectSum;
#endif
