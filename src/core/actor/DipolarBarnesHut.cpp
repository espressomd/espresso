#include "config.hpp"
#include "grid.hpp"
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

//std::cout << "Trace activate_dipolar_barnes_hut 1" << std::endl;
dipolarBarnesHut =new DipolarBarnesHut(espressoSystemInterface);
forceActors.push_back(dipolarBarnesHut);
energyActors.push_back(dipolarBarnesHut);
//std::cout << "Trace activate_dipolar_barnes_hut 2" << std::endl;

coulomb.Dmethod = DIPOLAR_BH_GPU;
}

void deactivate_dipolar_barnes_hut()
{
if (dipolarBarnesHut)
{
  //std::cout << "Trace deactivate_dipolar_barnes_hut 1" << std::endl;
  forceActors.remove(dipolarBarnesHut);
  energyActors.remove(dipolarBarnesHut);
  delete(dipolarBarnesHut);
  dipolarBarnesHut=NULL;
  //std::cout << "Trace deactivate_dipolar_barnes_hut 2" << std::endl;

}
}

DipolarBarnesHut *dipolarBarnesHut=0;

#endif

