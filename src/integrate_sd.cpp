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

/** \file integrate_sd.cpp   Stokes dynamics integrator.
 *
 *  For more information about the integrator 
 *  see \ref integrate_sd.hpp "integrate_sd.hpp".
*/

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "utils.hpp"
#include "integrate_sd.hpp"
#include "integrate_sd_cuda_debug.cuh"
#include "integrate.hpp"
#include "reaction.hpp"
#include "electrokinetics.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "communication.hpp"
#include "grid.hpp" // required
#include "cells.hpp"
#include "verlet.hpp"
#include "rotation.hpp"
#include "ghosts.hpp"
#include "pressure.hpp"
#include "p3m.hpp"
#include "maggs.hpp"
#include "thermostat.hpp"
#include "initialize.hpp"
#include "forces.hpp"
#include "nsquare.hpp"
#include "domain_decomposition.hpp"
#include "layered.hpp"
#include "nemd.hpp"
#include "rattle.hpp"
#include "errorhandling.hpp"
#include "lattice.hpp"
#include "lb.hpp"
#include "virtual_sites.hpp"
#include "statistics_correlation.hpp"
#include "ghmc.hpp"


/************************************************
 * DEFINES
 ************************************************/

/** Tag for communication in verlet fix: propagate_positions()  */
#define REQ_INT_VERLET   400

/*******************  variables  *******************/

//double time_step        = -1.0;

//double sim_time         = 0.0;

double sd_viscosity=-1;
double sd_radius=-1;

#ifdef ADDITIONAL_CHECKS
double db_max_force = 0.0, db_max_vel = 0.0;
int    db_maxf_id   = 0,   db_maxv_id = 0;
#endif

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Rescale all particle forces with <s>\f[ 0.5 \Delta t^2 \f]</s> some factor. */
//void rescale_forces_sd();
/** Propagate the positions. Integration step 2 of the Velocity Verletintegrator:<br>
    \f[ p(t+\Delta t) = p(t) + \Delta t  \mu v(t) \f] */
void propagate_pos_sd();

#ifdef CUDA
void propagate_pos_sd_cuda(double * box_l, int N,double * pos_h, double * force_h, double * velo_h);
#endif

int sd_get_particle_num();
/*@}*/

/************************************************************/


/************************************************************/

void integrate_sd(int n_steps)
{
  int i;

  /* Prepare the Integrator */
  on_integration_start();

  /* if any method vetoes (P3M not initialized), immediately bail out */
  if (check_runtime_errors())
    return;

  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv: integrating %d steps (recalc_forces=%d)\n",
		      this_node, n_steps, recalc_forces));
   
  /* Integration Step:
     Calculate forces f(t) as function of positions p(t) ( and velocities v(t) ) */
  //if (recalc_forces) { 
  //thermo_heat_up();
#ifdef LB
  transfer_momentum = 0;
  if (lattice_switch & LATTICE_LB && this_node == 0)
    fprintf (stderr, "Warning, no valid forces from previous integration step, means the LB fluid coupling is not included in the particle forces for this step.\n");
#endif
#ifdef LB_GPU
  transfer_momentum_gpu = 0;
  if (lattice_switch & LATTICE_LB_GPU && this_node == 0)
    fprintf (stderr, "Warning, no valid forces from previous integration step, means the GPU LB fluid coupling is not included in the particle forces for this step.\n");
#endif
//VIRTUAL_SITES pos (and vel for DPD) update for security reason !!!
#ifdef VIRTUAL_SITES
  update_mol_vel_pos();
  ghost_communicator(&cell_structure.update_ghost_pos_comm);
  if (check_runtime_errors()) return;
#endif
#ifdef COLLISION_DETECTION
  prepare_collision_queue();
#endif
  force_calc();
   
  //VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
  ghost_communicator(&cell_structure.collect_ghost_force_comm);
  init_forces_ghosts();
  distribute_mol_force();
  if (check_runtime_errors()) return;
#endif

  ghost_communicator(&cell_structure.collect_ghost_force_comm);

#ifdef ROTATION
  convert_initial_torques();
#endif

  //thermo_cool_down();

  /* Communication Step: ghost forces */


  /*apply trap forces to trapped molecules*/
#ifdef MOLFORCES
  // prob. works only with harmonic bounds
  calc_and_apply_mol_constraints();
#endif

  /* should be pretty late, since it needs to zero out the total force */
#ifdef COMFIXED
  calc_comfixed();
#endif

  //rescale_forces();
    
#ifdef COLLISION_DETECTION
  //should not be neccessery, as integrator avoids collision
  handle_collisions();
#endif
  // end of force calculation

#ifdef GHMC
  if(thermo_switch & THERMO_GHMC)
    ghmc_init();
#endif

  if (check_runtime_errors())
    return;

  n_verlet_updates = 0;

  /* Integration loop */
  for(i=0;i<n_steps;i++) {
    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,i));

#ifdef BOND_CONSTRAINT
    save_old_pos();
#endif

#ifdef GHMC
    if(thermo_switch & THERMO_GHMC) {
      if ((int) fmod(i,ghmc_nmd) == 0)
        ghmc_momentum_update();
    }
#endif
    if(thermo_switch & ~THERMO_SD){
      static bool warned_thermo_sd_other=false;
      if (!warned_thermo_sd_other){
	fprintf (stderr, "Warning, using another thermo than the one provided by StokesDynamics breaks (most likely) StokesDynamics.\n");
	warned_thermo_sd_other=true;
      }
    }

    /* Integration Steps: Update the Positions
       p_i(t + dt)   = p_i(t) + dt * mu_{ij} * f_j(t)
    */
    propagate_pos_sd(); // we dont have velocities

#ifdef BOND_CONSTRAINT
    const bool bond_constraint_with_sd_warned=false;
    if (!bond_constraint_with_sd_warned){ // warn only once
      fprintf (stderr, "Warning, using BOND_CONSTRAINT with StokesDynamics might not work as expected!.\n");    
      bond_constraint_with_sd_warned=true;
    }
    /**Correct those particle positions that participate in a rigid/constrained bond */
    cells_update_ghosts();

    correct_pos_shake();
#endif

#ifdef ELECTROSTATICS
    if(coulomb.method == COULOMB_MAGGS) {
      maggs_propagate_B_field(0.5*time_step); 
    }
#endif

#ifdef NPT
    if (check_runtime_errors())
      break;
#endif

    cells_update_ghosts();

    //VIRTUAL_SITES update pos and vel (for DPD)
#ifdef VIRTUAL_SITES
    update_mol_vel_pos();
    ghost_communicator(&cell_structure.update_ghost_pos_comm);

    if (check_runtime_errors()) break;
#if  defined(VIRTUAL_SITES_RELATIVE) && defined(LB) 
    // This is on a workaround stage: 
    // When using virtual sites relative and LB at the same time, it is necessary 
    // to reassemble the cell lists after all position updates, also of virtual
    // particles. 
    if (cell_structure.type == CELL_STRUCTURE_DOMDEC && (!dd.use_vList) ) 
      cells_update_ghosts();
#endif

#endif

    /* Integration Step: Step 3 of Velocity Verlet scheme:
       Calculate f(t+dt) as function of positions p(t+dt) ( and velocities v(t+0.5*dt) ) */
#ifdef LB
    transfer_momentum = 1;
#endif
#ifdef LB_GPU
    transfer_momentum_gpu = 1;
#endif

#ifdef COLLISION_DETECTION
    prepare_collision_queue();
#endif

    force_calc();

    //VIRTUAL_SITES distribute forces
#ifdef VIRTUAL_SITES
    ghost_communicator(&cell_structure.collect_ghost_force_comm);
    init_forces_ghosts();
    distribute_mol_force();
    if (check_runtime_errors()) break;
#endif

#ifdef CATALYTIC_REACTIONS
  integrate_reaction();
#endif

    /* Communication step: ghost forces */
    ghost_communicator(&cell_structure.collect_ghost_force_comm);

    /*apply trap forces to trapped molecules*/
#ifdef MOLFORCES         
    calc_and_apply_mol_constraints();
#endif

    /* should be pretty late, since it needs to zero out the total force */
#ifdef COMFIXED
    calc_comfixed();
#endif

    if (check_runtime_errors())
      break;

    
#ifdef LB
    if (lattice_switch & LATTICE_LB)
      lattice_boltzmann_update();
      
    if (check_runtime_errors())
      break;
#endif

#ifdef LB_GPU
    if(this_node == 0){
#ifdef ELECTROKINETICS
      if (ek_initialized) {
        ek_integrate();
      }
      else {
#endif
        if (lattice_switch & LATTICE_LB_GPU)
          lattice_boltzmann_update_gpu();
#ifdef ELECTROKINETICS
      }
#endif
    }
#endif //LB_GPU

#ifdef BOND_CONSTRAINT
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    correct_vel_shake();
#endif

#ifdef ROTATION
    convert_torques_propagate_omega();
#endif

    //VIRTUAL_SITES update vel
#ifdef VIRTUAL_SITES
    ghost_communicator(&cell_structure.update_ghost_pos_comm);
    update_mol_vel();
    if (check_runtime_errors()) break;
#endif

#ifdef ELECTROSTATICS
    if(coulomb.method == COULOMB_MAGGS) {
      maggs_propagate_B_field(0.5*time_step); 
    }
#endif

#ifdef COLLISION_DETECTION
    handle_collisions();
#endif

#ifdef NPT
    if((this_node==0) && (integ_switch == INTEG_METHOD_NPT_ISO))
      nptiso.p_inst_av += nptiso.p_inst;
#endif

#ifdef GHMC
    if(thermo_switch & THERMO_GHMC) {
      if ((int) fmod(i,ghmc_nmd) == ghmc_nmd-1)
        ghmc_mc();
    }
#endif

    /* Propagate time: t = t+dt */
    sim_time += time_step;
  }

  /* verlet list statistics */
  if(n_verlet_updates>0) verlet_reuse = n_steps/(double) n_verlet_updates;
  else verlet_reuse = 0;

#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    nptiso.invalidate_p_vel = 0;
    MPI_Bcast(&nptiso.p_inst, 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&nptiso.p_diff, 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast(&nptiso.volume, 1, MPI_DOUBLE, 0, comm_cart);
    if(this_node==0) nptiso.p_inst_av /= 1.0*n_steps;
    MPI_Bcast(&nptiso.p_inst_av, 1, MPI_DOUBLE, 0, comm_cart);
  }
#endif

#ifdef GHMC
  if(thermo_switch & THERMO_GHMC)
    ghmc_close();
#endif

}


/* Privat functions */
/************************************************************/


void propagate_pos_sd()
{
  /* Verlet list criterion */
  double skin2 = SQR(0.5 * skin);

  INTEG_TRACE(fprintf(stderr,"%d: propagate_pos:\n",this_node));
  Cell *cell;
  Particle *p;
  int c, i, np;
  //get total number of particles
  int N=sd_get_particle_num();
  // gather all the data for mobility calculation
  double * pos=NULL;
  pos=(double *)malloc(DIM*N*sizeof(double));
  assert(pos!=NULL);
  double * force=NULL;
  force=(double *)malloc(DIM*N*sizeof(double));
  assert(force!=NULL);
  double * velocity=NULL;
  bool old_velocity_reuse;
  { static int last_particle_number=0;
    if (N != last_particle_number){
      last_particle_number=N;
      velocity=(double *)calloc(DIM*N*sizeof(double),1);
      old_velocity_reuse=false;
    } else {
      old_velocity_reuse=true;
      velocity=(double *)malloc(DIM*N*sizeof(double));
    }}
  assert(velocity!=NULL);
#ifdef EXTERNAL_FORCES
  const int COORD_ALL=COORD_FIXED(0)&COORD_FIXED(1)&COORD_FIXED(2);
#endif
  int j=0; // total particle counter
  for (c = 0; c < local_cells.n; c++){
    cell = local_cells.cell[c];
    p    = cell->part;
    np   = cell->n;
    for (i = 0; i < np; i++) { // only count nonVirtual Particles
#ifdef EXTERNAL_FORCES
      if (p[i].l.ext_flag & COORD_ALL)
	{
	  fprintf (stderr, "Warning: Fixing particle in StokesDynamics this way with EXTERNAL_FORCES is not possible (and will be ignored). Please try to bind them e.g. harmonicaly.\n");
	}
#endif
#ifdef  VIRTUAL_SITES
      if (!ifParticleIsVirtual(&p[i]))
#endif
      {
        memcpy(&pos[3*j], p[i].r.p, 3*sizeof(double));
        memcpy(&force[3*j], p[i].f.f, 3*sizeof(double));
	if (old_velocity_reuse){
	  memcpy(velocity+3*j, p[i].m.v, 3*sizeof(double));
	}
	j++;
      }
    }
  }
  // cuda part
#ifdef CUDA
  //void propagate_pos_sd_cuda(double * box_l_h, int N,double * pos_h, double * force_h, double * velo_h){
  propagate_pos_sd_cuda(box_l,N,pos,force, velocity);
#else
  fprintf(stderr, "Warning - CUDA is currently required for SD");
#endif
  

#ifdef NEMD
  /* change momentum of each particle in top and bottom slab */
  fprintf (stderr, "Warning: NEMD is in SD not supported.\n");
#endif
  j=0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p    = cell->part;
    np   = cell->n;
    for (i = 0; i < np; i++) {
#ifdef VIRTUAL_SITES
      if (ifParticleIsVirtual(&p[i])) continue;
#endif
      // write back of position and velocity data
      memcpy(p[i].r.p, &pos[DIM*j], 3*sizeof(double));
      for (int d=0;d<DIM;d++){
	assert(!isnan(pos[DIM*i+d]));
      }
      memcpy(p[i].m.v, &velocity[DIM*j], 3*sizeof(double));
      j++;
      /* Verlet criterion check */
      if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) resort_particles = 1;
      
    }
  }
  /*std::cout << pos[0] << " " << pos[1] << " " << pos[2]
	    << pos[3] << " " << pos[4] << " " << pos[5]
	    << pos[6] << " " << pos[7] << " " << pos[8] << " "<< force << "\n";*/
  free(pos);
  free(force);
  announce_resort_particles();
  
}


int sd_set_particles_apart(){
  {
    double V=box_l[0];
    V*=box_l[1];
    V*=box_l[2];
    int N=sd_get_particle_num();
    double Vpart=4.*M_PI/3.*sd_radius*SQR(sd_radius)*N;
    double phi = Vpart/V;
    fprintf(stderr,"info: sd_set_particles_apart: Volume fraction of particles is %f.\n",phi);
    if (phi > M_PI/3./sqrt(2)){
      return -5;
    }
    if (phi > 0.7){
      return -4;
    }
    if (phi > 0.5){
      fprintf(stderr,"warning: sd_set_particles_apart: Volume fraction of particles is %f.\n\
It could be difficult to set particles apart!\n",phi);
    }
  }
  bool outer=false;
  do {
    outer=false;
    for (int c = 0; c < local_cells.n; c++){
      Cell * cell  = local_cells.cell[c];
      Particle * p = cell->part;
      int np   = cell->n;
      bool inner=false;
      do {
	inner=false;
	for (int i = 0; i < np; i++) { // only count nonVirtual Particles
	  for (int j = i+1; j <np; j++){
	    //position: p[i].r.p
	    double dr2=0;
	    double dr[3];
	    for (int d=0; d<3;d++){
	      dr[d]=p[i].r.p[d]-p[j].r.p[d];
	      dr[d]-=box_l[d]*rint(dr[d]/box_l[d]);
	      dr2+=SQR(dr[d]);
	    }
	    if (dr2 <= SQR(2*sd_radius*(1+1e-6))){
	      double drn = sqrt(dr2);
	      // push to a distance of 1e-6;
	      double fac=(sd_radius*(1+1.2e-6)-drn/2)/drn;
	      for (int d=0; d<3;d++){
		assert(!isnan(dr[d]*fac));
		p[i].r.p[d]+=fac*dr[d];
		p[j].r.p[d]-=fac*dr[d];
	      }
	      dr2=0;
	      for (int d=0; d<3;d++){
		dr[d]=p[i].r.p[d]-p[j].r.p[d];
		dr[d]-=box_l[d]*rint(dr[d]/box_l[d]);
		dr2+=SQR(dr[d]);
	      }
	      assert(dr2 > SQR(2*sd_radius*(1+1e-6)));
	      assert(dr2 < SQR(2*sd_radius*(1+2e-6)));
	      inner=true;
	    }
	  }
	}
      } while (inner);
      for (int cj =0; cj < local_cells.n;cj++){
	Cell * cellJ  = local_cells.cell[cj];
	Particle * pj = cellJ->part;
	int npj       = cellJ->n;
	for (int i = 0; i < np; i++) {
	  for (int j = (c==cj?i+1:0); j <npj; j++){
	    //position: p[i].r.p
	    double dr2=0;
	    double dr[3];
	    for (int d=0; d<3;d++){
	      dr[d]=p[i].r.p[d]-pj[j].r.p[d];
	      dr[d]-=box_l[d]*rint(dr[d]/box_l[d]);
	      dr2+=SQR(dr[d]);
	    }
	    if (dr2 <= SQR(2*sd_radius*(1+1e-6))){
	      double drn = sqrt(dr2);
	      // push to a distance of 1e-6;
	      double fac=(sd_radius*(1+1.2e-6)-drn/2)/drn;
	      assert(!isnan(fac));
	      assert (fac > 0);
	      for (int d=0; d<3;d++){
		p[i].r.p[d] +=fac*dr[d];
		pj[j].r.p[d]-=fac*dr[d];
	      }
	      outer=true;
	    }
	    for (int d=0;d<3;d++){
	      assert(!isnan( pj[j].r.p[d]));
	      assert(!isnan( p[i].r.p[d]));
	    }
	  }
	}
      }
    }
    fprintf(stderr,"+");
  } while (outer);
  fprintf(stderr,"set_apart suceeded ");
  return 0;
}


int sd_get_particle_num(){
  int N=0;
  for (int c = 0; c < local_cells.n; c++){
#ifdef VIRTUAL_SITES
    cell = local_cells.cell[c];
    p    = cell->part;
    np   = cell->n;
    for (i = 0; i < np; i++) { // only count nonVirtual Particles
      if (!ifParticleIsVirtual(&p[i])) ++N;
    } 
#else
    N  += local_cells.cell[c]->n;
#endif
  }
  return N;
}
