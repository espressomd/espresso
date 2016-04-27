/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
0*/

/** \file integrate_sd.cpp   Stokesian dynamics integrator.
 *
 *  For more information about the integrator 
 *  see \ref integrate_sd.hpp "integrate_sd.hpp".
*/

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <assert.h>
#include "utils.hpp"
#include "integrate_sd.hpp"
#include "integrate.hpp"
#include "reaction.hpp"
#include "electrokinetics.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "communication.hpp"
#include "grid.hpp"
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

using std::ostringstream;

/************************************************
 * DEFINES
 ************************************************/

/** Tag for communication in verlet fix: propagate_positions()  */
#define REQ_INT_VERLET   400

/*******************  variables  *******************/

/// switch between double and float
#ifdef SD_USE_FLOAT
typedef float real;
#else
typedef double real;
#endif

/*******************  variables  *******************/

double sd_viscosity=1;
double sd_radius=-1;

int sd_seed[]={0,('E'+'S'+'P'+'R'+'e'+'s'+'S'+'o')};
int sd_random_state[]={0,0};
double sd_random_precision=1e-3;

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Propagate the positions. Integration step:<br>
    \f[ p(t+\Delta t) = p(t) + \Delta t  \mu f(t) \f] */
void propagate_pos_sd();

#ifdef CUDA
void propagate_pos_sd_cuda(real * box_l, int N,real * pos_h, real * force_h, real * velo_h);
#endif

void propagate_pos_bd(int N, real * pos, real * force, real * velocity);

int sd_get_particle_num();
/*@}*/

void integrator_sanity_checks_sd()
{
  if ( time_step < 0.0 ) {
      ostringstream msg;
      msg <<"time_step not set";
      runtimeError(msg);
  }
  if ( temperature < 0.0 ) {
      ostringstream msg;
      msg <<"thermostat not initialized";
      runtimeError(msg);
  }
  if (sd_radius < 0) {
      ostringstream msg;
      msg <<"Stokesian Dynamics Hydrodynamic particle radius not initialized";
      runtimeError(msg);
  }
  if (sd_viscosity < 0) {
      ostringstream msg;
      msg <<"Stokesian Dynamics fluid viscosity not initialized";
      runtimeError(msg);
  }
}

/************************************************************/


/************************************************************/

#if defined(SD) || defined(BD)

void integrate_sd(int n_steps)
{
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
  for(int step=0;step<n_steps;step++) {
    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,step));
    //sd_set_particles_apart();
#ifdef BOND_CONSTRAINT
    save_old_pos();
#endif

#ifdef GHMC
    if(thermo_switch & THERMO_GHMC) {
      if (step % ghmc_nmd == 0)
        ghmc_momentum_update();
    }
#endif
    if(thermo_switch & ~(THERMO_SD|THERMO_BD) ){
      static bool warned_thermo_sd_other=false;
      if (!warned_thermo_sd_other){
	fprintf (stderr, "Warning, using another thermo than the one provided by StokesDynamics breaks (most likely) StokesDynamics.\n");
	warned_thermo_sd_other=true;
      }
    }
    if (thermo_switch & THERMO_SD && thermo_switch &THERMO_BD) {
      fprintf (stderr, "Warning: cannot use BD and SD. Disabeling BD!\n");
      thermo_switch &= ~THERMO_BD;
    }

    /* Integration Step: Step 3 of Velocity Verlet scheme:
       Calculate f(t) as function of positions p(t) ( and ``velocities'' v(t) ) */

#ifdef LB
    transfer_momentum = 1;
#endif
#ifdef LB_GPU
    transfer_momentum_gpu = 1;
#endif

    force_calc();

#ifdef CATALYTIC_REACTIONS
    integrate_reaction();
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

#ifdef ELECTROSTATICS
    if(coulomb.method == COULOMB_MAGGS) {
      maggs_propagate_B_field(0.5*time_step); 
    }
#endif

#ifdef NPT
    if((this_node==0) && (integ_switch == INTEG_METHOD_NPT_ISO))
      nptiso.p_inst_av += nptiso.p_inst;
#endif

#ifdef GHMC
    if(thermo_switch & THERMO_GHMC) {
      if (step % ghmc_nmd == ghmc_nmd-1)
        ghmc_mc();
    }
#endif



    /** Integration Steps: Update the Positions
      \[ p_i(t + dt)   = p_i(t) + dt * \mu_{ij} * f_j(t) + dt * \mu_{ij} * f^B_j \]
    */
    propagate_pos_sd(); // we dont have velocities

#ifdef BOND_CONSTRAINT
    static bool bond_constraint_with_sd_warned=false;
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

//#endif /* SD */

/* Privat functions */
/************************************************************/

//#ifdef SD

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
  real * pos=NULL;
  pos=(real *)Utils::malloc(DIM*N*sizeof(double));
  assert(pos!=NULL);
  real * force=NULL;
  force=(real *)Utils::malloc(DIM*N*sizeof(double));
  assert(force!=NULL);
  real * velocity=NULL;
  velocity=(real *)Utils::malloc(DIM*N*sizeof(real));
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
      if (p[i].p.ext_flag & COORD_ALL)
	{
	  fprintf (stderr, "Warning: Fixing particle in StokesDynamics this way with EXTERNAL_FORCES is not possible (and will be ignored). Please try to bind them e.g. harmonicaly.\n");
	}
#endif
#ifdef  VIRTUAL_SITES
      if (!ifParticleIsVirtual(&p[i]))
#endif
      {
#ifdef SD_USE_FLOAT
	for (int d=0;d<3;d++){
	  pos[3*j+d]        = p[i].r.p[d];
	  pos[3*j+d]        -=rint(pos[3*j+d]/box_l[d])*box_l[d];
	  force[3*j+d]      = p[i].f.f[d];
	}
#else
        memmove(&pos[3*j], p[i].r.p, 3*sizeof(double));
        memmove(&force[3*j], p[i].f.f, 3*sizeof(double));
	for (int d=0;d<3;d++){
	  pos[3*j+d]        -=rint(pos[3*j+d]/box_l[d])*box_l[d];
	}
#endif
	j++;
      }
    }
  }
  if (!(thermo_switch & THERMO_SD) && thermo_switch & THERMO_BD){
    propagate_pos_bd(N,pos,force, velocity);
  } else {
    // cuda part
#ifdef CUDA
    //void propagate_pos_sd_cuda(double * box_l_h, int N,double * pos_h, double * force_h, double * velo_h){
#ifdef SD_USE_FLOAT
    real box_size[3];
    for (int d=0;d<3;d++){
      box_size[d]=box_l[d];
    }
#else
    real * box_size = box_l;
#endif
    if (!(thermo_switch & THERMO_SD)){
      temperature*=-1;
    }
  propagate_pos_sd_cuda(box_size,N,pos,force, velocity);
  if (!(thermo_switch & THERMO_SD)){
    temperature*=-1;
  }
#else
  fprintf(stderr, "Warning - CUDA is currently required for SD\n");
  fprintf(stderr, "So i am just sitting here and copying stupidly stuff :'(\n");
#endif
}
  

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
#ifdef SD_USE_FLOAT
      for (int d=0;d<3;d++){
	p[i].r.p[d] = pos[3*j+d]+box_l[d]*rint(p[i].r.p[d]/box_l[d]);
	p[i].m.v[d] = velocity[3*j+d];
	//p[i].f.f[d] *= (0.5*time_step*time_step)/(*part).p.mass;
      }
#else
      for (int d=0;d<3;d++){
	p[i].r.p[d] = pos[3*j+d]+box_l[d]*rint(p[i].r.p[d]/box_l[d]);
      }
      memmove(p[i].m.v, &velocity[DIM*j], 3*sizeof(double));
#endif
      // somehow this does not effect anything, although it is called ...
      for (int d=0;d<3;d++){
	p[i].f.f[d] *= (0.5*time_step*time_step)/(*p).p.mass;
      }
      for (int d=0;d<DIM;d++){
        assert(!isnan(pos[DIM*i+d]));
      }
      j++;
      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PPOS p = (%.3f,%.3f,%.3f)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2]));
   

#ifdef ROTATION
      propagate_omega_quat_particle(&p[i]);
#endif

      /* Verlet criterion check */
      if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) resort_particles = 1;
    }
  }
  free(pos);
  free(force);
  free(velocity);
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
    //fprintf(stderr,"info: sd_set_particles_apart: Volume fraction of particles is %f.\n",phi);
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
  bool moved=false;
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
	    if (dr2 <= SQR(2*sd_radius*(1+1e-5))){
	      double drn = sqrt(dr2);
	      // push to a distance of 1e-5;
	      double fac=(sd_radius*(1+1.2e-5)-drn/2)/drn;
	      assert(!isnan(fac));
	      //printf("%d %d\t\t",i,j);
	      for (int d=0; d<3;d++){
		if(isnan(dr[d]*fac)){
		  fprintf(stderr, "%4d %4d %6e %6e ",i,j,dr[d],fac);
		}
		assert(!isnan(dr[d]));
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
	      assert(dr2 > SQR(2*sd_radius*(1+1e-5)));
	      assert(dr2 < SQR(2*sd_radius*(1+2e-5)));
	      inner=true;
	      moved=true;
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
	    if (dr2 <= SQR(2*sd_radius*(1+1e-5))){
	      double drn = sqrt(dr2);
	      // push to a distance of 1e-5;
	      double fac=(sd_radius*(1+1.2e-5)-drn/2)/drn;
	      assert(!isnan(fac));
	      assert (fac > 0);
	      for (int d=0; d<3;d++){
		p[i].r.p[d] +=fac*dr[d];
		pj[j].r.p[d]-=fac*dr[d];
	      }
	      outer=true;
	      moved=true;
	    }
	    for (int d=0;d<3;d++){
	      assert(!isnan( pj[j].r.p[d]));
	      assert(!isnan( p[i].r.p[d]));
	    }
	  }
	}
      }
    }
    if (moved) fprintf(stderr,"+");
  } while (outer);
  //fprintf(stderr,"set_apart suceeded ");
  for (int c = 0; c < local_cells.n; c++){
    Cell * cell  = local_cells.cell[c];
    Particle * p = cell->part;
    int np   = cell->n;
    for (int cj =0; cj < local_cells.n;cj++){
      Cell * cellJ  = local_cells.cell[cj];
      Particle * pj = cellJ->part;
      int npj       = cellJ->n;
      for (int i = 0; i < np; i++) {
	for (int j = (c==cj?i+1:0); j <npj; j++){
	  double dr2=0;
	  double dr[3];
	  for (int d=0; d<3;d++){
	    dr[d]=p[i].r.p[d]-pj[j].r.p[d];
	    dr[d]-=box_l[d]*rint(dr[d]/box_l[d]);
	    dr2+=SQR(dr[d]);
	  }
	  assert((dr2 > SQR(2*sd_radius*(1+1e-5))));
	}
      }
    }
  }
  return 0;
}

int sd_get_particle_num(){
  int N=0;
  for (int c = 0; c < local_cells.n; c++){
#ifdef VIRTUAL_SITES
    Cell * cell = local_cells.cell[c];
    Particle * p    = cell->part;
    int np   = cell->n;
    for (int i = 0; i < np; i++) { // only count nonVirtual Particles
      if (!ifParticleIsVirtual(&p[i])) ++N;
    } 
#else
    N  += local_cells.cell[c]->n;
#endif
  }
  return N;
}

void propagate_pos_bd(int N, real * pos, real * force, real * velocity){
  real self_mobility=1/(6*M_PI*sd_viscosity*sd_radius);
  real scal_f = time_step*self_mobility;
  real scal_r=sqrt(2*temperature*time_step*self_mobility);
  for (int i=0;i<3*N;i++){
    pos[i]+=scal_f * force[i]
          + scal_r * gaussian_random();
  }
}

#endif /* SD */
