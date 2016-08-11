#include "config.hpp"
#include "grid.hpp"
#include "DipolarBarnesHut_cuda.hpp"
#include "DipolarBarnesHut.hpp"
#include "../forces.hpp"
#include "EspressoSystemInterface.hpp"
#include "forces.hpp"
#include "energy.hpp"

#ifdef BARNES_HUT

void activate_dipolar_barnes_hut()
{
if (dipolarBarnesHut)
  free(dipolarBarnesHut);

dipolarBarnesHut =new DipolarBarnesHut(espressoSystemInterface);
forceActors.push_back(dipolarBarnesHut);
energyActors.push_back(dipolarBarnesHut);

coulomb.Dmethod = DIPOLAR_BH_GPU;
}

void deactivate_dipolar_barnes_hut()
{
if (dipolarBarnesHut)
{
  forceActors.remove(dipolarBarnesHut);
  energyActors.remove(dipolarBarnesHut);
  delete(dipolarBarnesHut);
  dipolarBarnesHut=NULL;

}
}

DipolarBarnesHut *dipolarBarnesHut=0;

#endif

