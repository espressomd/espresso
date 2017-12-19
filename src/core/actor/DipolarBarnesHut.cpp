#include "config.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "DipolarBarnesHut.hpp"
#include "../forces.hpp"
#include "EspressoSystemInterface.hpp"
#include "forces.hpp"
#include "energy.hpp"

#ifdef DIPOLAR_BARNES_HUT

void activate_dipolar_barnes_hut(float epssq, float itolsq)
{
    if (dipolarBarnesHut)
        free(dipolarBarnesHut);

    // also necessary on 1 CPU or GPU, does more than just broadcasting
    mpi_bcast_coulomb_params();
    dipolarBarnesHut =new DipolarBarnesHut(espressoSystemInterface, epssq, itolsq);
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

