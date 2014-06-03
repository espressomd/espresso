
/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#include <mpi.h>
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

#include "EspressoSystemInterface.hpp"

PotentialList potentials;

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
  // Communication step: distribute ghost positions
  cells_update_ghosts();

  // VIRTUAL_SITES pos (and vel for DPD) update for security reason !!!
#ifdef VIRTUAL_SITES
  update_mol_vel_pos();
  ghost_communicator(&cell_structure.update_ghost_pos_comm);
#endif

#if defined(VIRTUAL_SITES_RELATIVE) && defined(LB) 
  // This is on a workaround stage: 
  // When using virtual sites relative and LB at the same time, it is necessary 
  // to reassemble the cell lists after all position updates, also of virtual
  // particles. 
  if ((lattice_switch & LATTICE_LB) && cell_structure.type == CELL_STRUCTURE_DOMDEC && (!dd.use_vList) ) 
    cells_update_ghosts();
#endif
  
#ifdef COLLISION_DETECTION
  prepare_collision_queue();
#endif

  espressoSystemInterface.update();

  // Compute the forces from the force objects
  for (PotentialList::iterator potential= potentials.begin();
		  potential != potentials.end(); ++potential)
	  (*potential)->computeForces(espressoSystemInterface);

#ifdef LB_GPU
#ifdef SHANCHEN
  if (lattice_switch & LATTICE_LB_GPU && this_node == 0) lattice_boltzmann_calc_shanchen_gpu();
#endif // SHANCHEN

  // transfer_momentum_gpu check makes sure the LB fluid doesn't get updated on integrate 0
  // this_node==0 makes sure it is the master node where the gpu exists
  if (lattice_switch & LATTICE_LB_GPU && transfer_momentum_gpu && (this_node == 0) ) lb_calc_particle_lattice_ia_gpu();
#endif // LB_GPU

#ifdef ELECTROSTATICS
  if (iccp3m_initialized && iccp3m_cfg.set_flag)
    iccp3m_iteration();
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

#ifdef VOLUME_FORCE
	double volume=0.;
	
	for (int i=0;i< MAX_OBJECTS_IN_FLUID;i++){
		calc_volume(&volume,i);
		if (volume<1e-100) break;
		add_volume_force(volume,i);	
	}
#endif	

#ifdef AREA_FORCE_GLOBAL
	double area=0.;

	for (int i=0;i< MAX_OBJECTS_IN_FLUID;i++){
		calc_area_global(&area,i);
		if (area<1e-100) break;
		add_area_global_force(area,i);
	}
#endif	

  calc_long_range_forces();

#ifdef LB
  if (lattice_switch & LATTICE_LB) calc_particle_lattice_ia() ;
#endif

#ifdef COMFORCE
  calc_comforce();
#endif

#ifdef METADYNAMICS
  /* Metadynamics main function */
  meta_perform();
#endif

#ifdef CUDA
  copy_forces_from_GPU();
#endif

  // VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
  ghost_communicator(&cell_structure.collect_ghost_force_comm);
  init_forces_ghosts();
  distribute_mol_force();
#endif

  // Communication Step: ghost forces
  ghost_communicator(&cell_structure.collect_ghost_force_comm);

  // apply trap forces to trapped molecules
#ifdef MOLFORCES         
  calc_and_apply_mol_constraints();
#endif

  // should be pretty late, since it needs to zero out the total force
#ifdef COMFIXED
  calc_comfixed();
#endif

  // mark that forces are now up-to-date
  recalc_forces = 0;

#ifdef COLLISION_DETECTION
  handle_collisions();
#endif
}

/************************************************************/

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
#ifdef CUDA
    case COULOMB_P3M_GPU:
      if (this_node == 0) p3m_gpu_add_farfield_force();
  #ifdef NPT
      printf("NPT can not be used in conjunction with the GPU P3M\n"); //TODO fix this?
      exit(1); //TODO ALTERNATIVELY CHECK IF BAROSTAT IS ACTUALLY ON
  #endif
      break;
#endif
    case COULOMB_P3M:
      p3m_charge_assign();
  #ifdef NPT
      if(integ_switch == INTEG_METHOD_NPT_ISO)
        nptiso.p_vir[0] += p3m_calc_kspace_forces(1,1);
      else
  #endif
        p3m_calc_kspace_forces(1,0);
      break;
  #endif
    case COULOMB_MAGGS:
      maggs_calc_forces();
      break;
    case COULOMB_MMM2D:
      MMM2D_add_far_force();
      MMM2D_dielectric_layers_force_contribution();
    }
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
  }
#endif  /*ifdef DIPOLES */
}

/************************************************************/

/** initialize the forces for a real particle */
inline void init_local_particle_force(Particle *part)
{
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

    #ifdef EXTERNAL_FORCES
      if(part->l.ext_flag & PARTICLE_EXT_TORQUE) {
        part->f.torque[0] += part->l.ext_torque[0];
        part->f.torque[1] += part->l.ext_torque[1];
        part->f.torque[2] += part->l.ext_torque[2];
      }
    #endif
    
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

/** initialize the forces for a ghost particle */
inline void init_ghost_force(Particle *part)
{
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


