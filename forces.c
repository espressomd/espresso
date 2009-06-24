// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
/** \file forces.c Force calculation.
 *
 *  For more information see \ref forces.h "forces.h".
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "thermostat.h"
#include "pressure.h"
#include "communication.h"
#include "ghosts.h" 
#include "verlet.h"
#include "grid.h"
#include "cells.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "rotation.h"
#include "forces.h"
#include "elc.h"
#include "lattice.h"
#include "lb.h"
#include "nsquare.h"
#include "layered.h"
#include "domain_decomposition.h"
#include "magnetic_non_p3m__methods.h"
#include "mdlc_correction.h"
#include "virtual_sites.h"

/************************************************************/
/* local prototypes                                         */
/************************************************************/

/** Calculate long range forces (P3M, MMM2d...). */
void calc_long_range_forces();

/** initialize real particle forces with thermostat forces and
    ghost particle forces with zero. */
void init_forces();

/************************************************************/

void force_calc()
{
  

#if defined(DIPOLES) && defined(ROTATION)
  convert_quat_to_dip_all();
#endif

  init_forces();
  
  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    layered_calculate_ia();
    break;
  case CELL_STRUCTURE_DOMDEC:
    if(dd.use_vList) {
      if (rebuild_verletlist)
	build_verlet_lists_and_calc_verlet_ia();
      else
	calculate_verlet_ia();
    }
    else
      calc_link_cell();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_calculate_ia();
    
  }
  
  calc_long_range_forces();

#ifdef LB
  if (lattice_switch & LATTICE_LB) calc_particle_lattice_ia() ;
#endif

#ifdef COMFORCE
  calc_comforce();
#endif

/* this must be the last force to be calculated (Mehmet)*/
#ifdef COMFIXED
  calc_comfixed();
#endif

}

/************************************************************/

void calc_long_range_forces()
{
#ifdef ELECTROSTATICS  
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
#ifdef ELP3M
  case COULOMB_ELC_P3M:
    if (elc_params.dielectric_contrast_on) {
      ELC_P3M_modify_p3m_sums_both();
      ELC_P3M_charge_assign_both();
      ELC_P3M_self_forces();
    }
    else
      P3M_charge_assign();

    P3M_calc_kspace_forces_for_charges(1,0);

    if (elc_params.dielectric_contrast_on)
      ELC_P3M_restore_p3m_sums();
 
    ELC_add_force(); 

    break;
  case COULOMB_P3M:
    P3M_charge_assign();
#ifdef NPT
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += P3M_calc_kspace_forces_for_charges(1,1);
    else
#endif
      P3M_calc_kspace_forces_for_charges(1,0);
    break;
#endif
  case COULOMB_EWALD:
#ifdef NPT
    if(integ_switch == INTEG_METHOD_NPT_ISO)
      nptiso.p_vir[0] += EWALD_calc_kspace_forces(1,1);
    else
#endif
      EWALD_calc_kspace_forces(1,0);
    break;
  case COULOMB_MAGGS:
    maggs_calc_e_forces();
    break;
  case COULOMB_MMM2D:
    MMM2D_add_far_force();
    MMM2D_dielectric_layers_force_contribution();
  }
#endif  /*ifdef ELECTROSTATICS */

#ifdef MAGNETOSTATICS  
  /* calculate k-space part of the magnetostatic interaction. */
  switch (coulomb.Dmethod) {
#ifdef ELP3M
#ifdef MDLC
  case DIPOLAR_MDLC_P3M:
     add_mdlc_force_corrections();
    //fall through 
#endif
  case DIPOLAR_P3M:
    P3M_dipole_assign();
#ifdef NPT
    if(integ_switch == INTEG_METHOD_NPT_ISO) {
      nptiso.p_vir[0] += P3M_calc_kspace_forces_for_dipoles(1,1);
      fprintf(stderr,"dipolar_P3M at this moment is added to p_vir[0]\n");    
    } else
#endif
      P3M_calc_kspace_forces_for_dipoles(1,0);

      break;
#endif
#ifdef DAWAANR
  case DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA: 
      dawaanr_calculations(1,0);
      break;
#endif
#ifdef MAGNETIC_DIPOLAR_DIRECT_SUM
#ifdef MDLC
  case DIPOLAR_MDLC_DS:
     add_mdlc_force_corrections();
    //fall through 
#endif
  case DIPOLAR_DS: 
        magnetic_dipolar_direct_sum_calculations(1,0);
      break;
#endif

  }
#endif  /*ifdef MAGNETOSTATICS */
}

/************************************************************/

/** initialize the forces for a real particle */
MDINLINE void init_local_particle_force(Particle *part)
{
  //  if (part->p.identity == 9)
  // {
  //  printf("particle 9 local, weight %f  ", part->p.adress_weight);
#ifdef ADRESS
  if (ifParticleIsVirtual(part)) {
    part->p.adress_weight=adress_wf_vector(part->r.p);
  }
#endif
  if ( thermo_switch & THERMO_LANGEVIN )
    friction_thermo_langevin(part);
  else {
    part->f.f[0] = 0;
    part->f.f[1] = 0;
    part->f.f[2] = 0;
  }

#ifdef EXTERNAL_FORCES   
  if(part->l.ext_flag & PARTICLE_EXT_FORCE) {
    part->f.f[0] += part->l.ext_force[0];
    part->f.f[1] += part->l.ext_force[1];
    part->f.f[2] += part->l.ext_force[2];
  }
#endif
  
#ifdef ROTATION
  {
    double scale;
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;
    
    /* and rescale quaternion, so it is exactly of unit length */	
    scale = sqrt( SQR(part->r.quat[0]) + SQR(part->r.quat[1]) +
		  SQR(part->r.quat[2]) + SQR(part->r.quat[3]));
    part->r.quat[0]/= scale;
    part->r.quat[1]/= scale;
    part->r.quat[2]/= scale;
    part->r.quat[3]/= scale;
  }
#endif

#ifdef ADRESS
  /** #ifdef THERMODYNAMIC_FORCE */
  if(ifParticleIsVirtual(part))
    if(part->p.adress_weight > 0 && part->p.adress_weight < 1)
      add_thermodynamic_force(part);
  /** #endif */  
#endif
}

/** initialize the forces for a ghost particle */
MDINLINE void init_ghost_force(Particle *part)
{
  //if (part->p.identity == 9)
  // {
  //  printf("particle 9 ghost, weight %f  ", part->p.adress_weight);
  #ifdef ADRESS
  if (ifParticleIsVirtual(part)) {
    part->p.adress_weight=adress_wf_vector(part->r.p);
  }
  #endif
  
  part->f.f[0] = 0;
  part->f.f[1] = 0;
  part->f.f[2] = 0;

#ifdef ROTATION
  {
    double scale;
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;

    /* and rescale quaternion, so it is exactly of unit length */	
    scale = sqrt( SQR(part->r.quat[0]) + SQR(part->r.quat[1]) +
		  SQR(part->r.quat[2]) + SQR(part->r.quat[3]));
    part->r.quat[0]/= scale;
    part->r.quat[1]/= scale;
    part->r.quat[2]/= scale;
    part->r.quat[3]/= scale;
  }
#endif
}

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


