#include "config.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
#include "VirtualSitesInertialessTracers.hpp"
#include "cells.hpp"
#include "errorhandling.hpp"
#include "lattice.hpp"
#include "virtual_sites/lb_inertialess_tracers.hpp"

void VirtualSitesInertialessTracers::after_force_calc() {
  // Now the forces are computed and need to go into the LB fluid
#ifdef LB
  if (lattice_switch & LATTICE_LB) {
    IBM_ForcesIntoFluid_CPU();
    return;
  }
#endif
#ifdef LB_GPU
  if (lattice_switch & LATTICE_LB_GPU) {
    IBM_ForcesIntoFluid_GPU(local_cells.particles());
    return;
  }
#endif
  runtimeErrorMsg() << "Inertialess Tracers: No LB method was active.";
}

void VirtualSitesInertialessTracers::after_lb_propagation() {
#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

  IBM_UpdateParticlePositions(local_cells.particles());
// We reset all since otherwise the halo nodes may not be reset
// NB: the normal Espresso reset is also done after applying the forces
//    if (lattice_switch & LATTICE_LB) IBM_ResetLBForces_CPU();
#ifdef LB_GPU
// if (lattice_switch & LATTICE_LB_GPU) IBM_ResetLBForces_GPU();
#endif

  // Ghost positions are now out-of-date
  // We should update.
  // Actually we seem to get the same results whether we do this here or not,
  // but it is safer to do it
  ghost_communicator(&cell_structure.update_ghost_pos_comm);

#endif // VS inertialess tracers
}
#endif
