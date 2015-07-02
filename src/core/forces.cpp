
/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file forces.cpp Force calculation.
 *
 *  For more information see \ref forces.hpp "forces.h".
*/
/*#include <mpi.h>
#include "forces_inline.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "utils.hpp"
#include "thermostat.hpp"
#include "pressure.hpp"
#include "communication.hpp"
#include "ghosts.hpp" 
#include "verlet.hpp"
#include "grid.hpp"
#include "cells.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "rotation.hpp"
#include "forces.hpp"
#include "elc.hpp"
#include "lattice.hpp"
#include "lb.hpp"
#include "nsquare.hpp"
#include "layered.hpp"
#include "domain_decomposition.hpp"
#include "magnetic_non_p3m_methods.hpp"
#include "mdlc_correction.hpp"
#include "virtual_sites.hpp"
#include "constraint.hpp"
#include "lbgpu.hpp"
#include "iccp3m.hpp"
#include "p3m_gpu.hpp"
#include "cuda_interface.hpp"

#include "EspressoSystemInterface.hpp"*/

#include "p3m_gpu.hpp"
#include "maggs.hpp"
#include "forces_inline.hpp"
#include "electrokinetics.hpp"
ActorList forceActors;

void init_forces()
{
  Cell *cell;
  Particle *p;
  int np, c, i;

  /* The force initialization depends on the used thermostat and the
     thermodynamic ensemble */

#ifdef NPT
  /* reset virial part of instantaneous pressure */
  if(integ_switch == INTEG_METHOD_NPT_ISO) 
    nptiso.p_vir[0] = nptiso.p_vir[1] = nptiso.p_vir[2] = 0.0;
#endif


  /* initialize forces with langevin thermostat forces
     or zero depending on the thermostat
     set torque to zero for all and rescale quaternions
  */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++)
      init_local_particle_force(&p[i]);
  }

  /* initialize ghost forces with zero
     set torque to zero for all and rescale quaternions
  */
  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++)
      init_ghost_force(&p[i]);
  }

#ifdef CONSTRAINTS
  init_constraint_forces();
#endif
}

void init_forces_ghosts()
{
  Cell *cell;
  Particle *p;
  int np, c, i;

  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++)
      init_ghost_force(&p[i]);
  }
}

// This function is no longer called from force_calc().
// The check was moved to rescale_fores() to avoid an additional iteration over all particles
void check_forces()
{
  Cell *cell;
  Particle *p;
  int np, c, i;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++) {
      check_particle_force(&p[i]);
    }
  }

  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++)
      check_particle_force(&p[i]);
  }
}

void calc_long_range_forces()
{
#ifdef ELECTROSTATICS  
	/* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
#ifdef P3M
  case COULOMB_ELC_P3M:
    if (elc_params.dielectric_contrast_on) {
      ELC_P3M_modify_p3m_sums_both();
      ELC_p3m_charge_assign_both();
      ELC_P3M_self_forces();
    }
    else
      p3m_charge_assign();
    
    p3m_calc_kspace_forces(1,0);
    
    if (elc_params.dielectric_contrast_on)
      ELC_P3M_restore_p3m_sums();
    
    ELC_add_force();
    
    break;
#endif
#ifdef CUDA
  case COULOMB_P3M_GPU:
    if (this_node == 0) {
      FORCE_TRACE(printf("Computing GPU P3M forces.\n"));
      p3m_gpu_add_farfield_force();
    }
    /* there is no NPT handling here as long as we cannot compute energies.
       This is checked in integrator_npt_sanity_checks() when integration starts. */
    break;
#endif
#ifdef P3M
  case COULOMB_P3M:
    FORCE_TRACE(printf("%d: Computing P3M forces.\n", this_node));
    p3m_charge_assign();
#ifdef NPT
    if (integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += p3m_calc_kspace_forces(1,1);
    else
#endif
      p3m_calc_kspace_forces(1, 0);
    break;
#endif
  case COULOMB_MAGGS:
    maggs_calc_forces();
    break;
  case COULOMB_MMM2D:
    MMM2D_add_far_force();
    MMM2D_dielectric_layers_force_contribution();
    break;
  default:
    break;
  }

/* If enabled, calculate electrostatics contribution from electrokinetics species. */ 
#ifdef EK_ELECTROSTATIC_COUPLING
  ek_calculate_electrostatic_coupling();
#endif

#endif  /*ifdef ELECTROSTATICS */
 
#ifdef DIPOLES  
  /* calculate k-space part of the magnetostatic interaction. */
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    add_mdlc_force_corrections();
    //fall through 
  case DIPOLAR_P3M:
    dp3m_dipole_assign();
#ifdef NPT
    if(integ_switch == INTEG_METHOD_NPT_ISO) {
      nptiso.p_vir[0] += dp3m_calc_kspace_forces(1,1);
      fprintf(stderr,"dipolar_P3M at this moment is added to p_vir[0]\n");    
    } else
#endif
      dp3m_calc_kspace_forces(1,0);
    
    break;
#endif
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA: 
    dawaanr_calculations(1,0);
    break;
  case DIPOLAR_MDLC_DS:
    add_mdlc_force_corrections();
    //fall through 
  case DIPOLAR_DS: 
    magnetic_dipolar_direct_sum_calculations(1,0);
    break;
  case DIPOLAR_DS_GPU: 
    // Do nothing. It's an actor
    break;
  case DIPOLAR_NONE:
      break;
  default:
      ostringstream msg;
      msg <<"unknown dipolar method";
      runtimeError(msg);
      break;
  }
#endif  /*ifdef DIPOLES */
}

void
calc_non_bonded_pair_force_from_partcfg(Particle *p1, Particle *p2, IA_parameters *ia_params,
                                        double d[3], double dist, double dist2,
                                        double force[3],
                                        double torque1[3], double torque2[3]) {
#ifdef MOL_CUT
   //You may want to put a correction factor and correction term for smoothing function else then theta
   if (checkIfParticlesInteractViaMolCut_partcfg(p1,p2,ia_params)==1)
#endif
   {
     calc_non_bonded_pair_force_parts(p1, p2, ia_params,
                                      d, dist, dist2, force, torque1, torque2);
   }
}

void
calc_non_bonded_pair_force_from_partcfg_simple(Particle *p1, Particle *p2,
                                               double d[3], double dist,
                                               double dist2, double force[3]){
   IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
   double torque1[3],torque2[3];
   calc_non_bonded_pair_force_from_partcfg(p1, p2, ia_params, d, dist, dist2,
                                           force, torque1, torque2);
}


