// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file integrate.c   Molecular dynamics integrator.
 *
 *  For more information about the integrator 
 *  see \ref integrate.h "integrate.h".
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "integrate.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "communication.h"
#include "grid.h"
#include "cells.h"
#include "verlet.h"
#include "rotation.h"
#include "ghosts.h"
#include "debug.h"
#include "pressure.h"
#include "p3m.h"
#include "utils.h"
#include "thermostat.h"
#include "initialize.h"
#include "forces.h"
#include "nsquare.h"
#include "domain_decomposition.h"
#include "layered.h"
#include "nemd.h"

/************************************************
 * DEFINES
 ************************************************/

/** Tag for communication in verlet fix: propagate_positions()  */
#define REQ_INT_VERLET   400

/*******************  variables  *******************/

int integ_switch;

double time_step = -1.0;
double sim_time  = 0.0;

double skin       = -1.0;
double skin2;
double max_range  = -1.0;
double max_range2 = -1.0;

int    resort_particles = 1;
int    recalc_forces = 1;

double verlet_reuse=0.0;

#ifdef ADDITIONAL_CHECKS
double db_max_force=0.0, db_max_vel=0.0;
int    db_maxf_id=0,db_maxv_id=0;
#endif

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Rescale all particle forces with \f[ 0.5 \Delta t^2 \f]. */
void rescale_forces();
/** Propagate the velocities and positions */
void propagate_vel_pos();
/** Rescale all particle forces with \f[ 0.5 \Delta t^2 \f] and propagate the velocities */
void rescale_forces_propagate_vel();

void force_and_velocity_check(Particle *p); 

void force_and_velocity_display();
 
void finalize_p_inst_npt();

/*@}*/

/************************************************************/

int invalidate_system(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  mpi_bcast_event(INVALIDATE_SYSTEM);
  return TCL_OK;
}

void local_invalidate_system()
{
  resort_particles = 1;
  invalidate_obs();
}

int integrate(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int  n_steps;
  
  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <step num> \"", (char *) NULL);
    return (TCL_ERROR);
  }

  /* translate argument */
  if (Tcl_GetInt(interp, argv[1], &n_steps) == TCL_ERROR)
    return (TCL_ERROR);

  if(n_steps < 0) {
    Tcl_AppendResult(interp, "illegal number of steps", (char *) NULL);
    return (TCL_ERROR);
  }

  /* assume velocity verlet integration with langevin thermostat */
  if (argc != 2) {
    Tcl_AppendResult(interp, "too many arguments:  should be \"",
		     argv[0], " <task> \"", (char *) NULL);
    return (TCL_ERROR);
  }

  mpi_integrate(n_steps);

  return (TCL_OK);
}

/************************************************************/

void integrate_vv_recalc_maxrange()
{
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_recalc_maxrange:\n",this_node));

  /* maximal interaction cutoff */
  calc_maximal_cutoff();
  max_range  = max_cut;
  /* at beginning be nice */
  if (skin > 0.0)
    max_range += skin;
  max_range2 = max_range * max_range;
}

/************************************************************/
void integrate_ensemble_init()
{
#ifdef NPT
  if (piston != 0.0) {
    /* prepare NpT-integration */
    inv_piston = 1/(1.0*piston);
    NpT_volume = box_l[0]*box_l[1]*box_l[2];
    if (recalc_forces) { 
      p_inst = 0.0;  
      p_vir  = 0.0; 
      p_vel  = 0.0; 
    }
  }
#endif
}

/************************************************************/

void integrate_vv(int n_steps)
{
  int i;
  int n_verlet_updates = 0;

  /* Prepare the Integrator */
  on_integration_start();

  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv: integrating %d steps (recalc_forces=%d)\n",
		      this_node, n_steps, recalc_forces));
   
  /* Integration Step: Preparation for first integration step:
     Calculate forces f(t) as function of positions p(t) ( and velocities v(t) ) */
  if (recalc_forces) {
    force_calc(); 
#ifdef ROTATION
    convert_initial_torques();
#endif

    /* Communication Step: ghost forces */
    ghost_communicator(&cell_structure.collect_ghost_force_comm);
    rescale_forces();
    recalc_forces = 0;
  }

  /* Integration loop */
  for(i=0;i<n_steps;i++) {

    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,i));

    /* Integration Steps: Step 1 and 2 of Velocity Verlet scheme:
       v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
       p(t + dt)   = p(t) + dt * v(t+0.5*dt)
       NOTE: Prefactors do not occur in formulas since we use 
       rescaled forces and velocities. */
    propagate_vel_pos();
#ifdef ROTATION
    propagate_omega_quat(); 
#endif

    /* Check verlet list criterion */
    if(rebuild_verletlist == 1 &&
       cell_structure.type == CELL_STRUCTURE_DOMDEC) {
      INTEG_TRACE(fprintf(stderr,"%d: Rebuild Verlet List\n",this_node));
      n_verlet_updates++;
      /* Communication step:  number of ghosts and ghost information */
      cells_resort_particles(CELL_NEIGHBOR_EXCHANGE);
    }
    else {
      /* Communication step: ghost information */
      ghost_communicator(&cell_structure.update_ghost_pos_comm);
    }

    /* Integration Step: Step 3 of Velocity Verlet scheme:
       Calculate f(t+dt) as function of positions p(t+dt) ( and velocities v(t+0.5*dt) ) */
    force_calc();

    /* Communication step: ghost forces */
    ghost_communicator(&cell_structure.collect_ghost_force_comm);

    /* Integration Step: Step 4 of Velocity Verlet scheme:
       v(t+dt) = v(t+0.5*dt) + 0.5*dt * f(t+dt) */
    rescale_forces_propagate_vel();
#ifdef ROTATION
    convert_torqes_propagate_omega();
#endif        

    /* Propagate time: t = t+dt */
    if(this_node==0) sim_time += time_step;
  }

  /* verlet list statistics */
  if(n_verlet_updates>0) verlet_reuse = n_steps/(double) n_verlet_updates;
  else verlet_reuse = 0;

#ifdef NPT
  if (piston != 0.0) {
    MPI_Bcast(&p_inst,     1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p_diff,     1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NpT_volume, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
#endif
}

/************************************************************/

void rescale_velocities(double scale) 
{
  Particle *p;
  int i, np, c;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      p[i].m.v[0] *= scale;
      p[i].m.v[1] *= scale;
      p[i].m.v[2] *= scale;
    }
  }
}

/* Callback functions */
/************************************************************/


int skin_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "skin must be positive.", (char *) NULL);
    return (TCL_ERROR);
  }
  skin = data;
  mpi_bcast_parameter(FIELD_SKIN);
  return (TCL_OK);
}

int time_step_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0.0) {
    Tcl_AppendResult(interp, "time step must be positive.", (char *) NULL);
    return (TCL_ERROR);
  }
  time_step = data;
  mpi_set_time_step();

  return (TCL_OK);
}

int time_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  sim_time = data;
  mpi_bcast_parameter(FIELD_SIMTIME);
  return (TCL_OK);
}


/* Privat functions */
/************************************************************/

void rescale_forces()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  double scale ;

  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces:\n",this_node));

  scale = 0.5 * time_step * time_step;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      p[i].f.f[0] *= scale;
      p[i].f.f[1] *= scale;
      p[i].f.f[2] *= scale;

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));

    }
  }
}

void rescale_forces_propagate_vel() 
{
  Cell *cell;
  Particle *p;
  int i, j, np, c;
  double scale;

#ifdef NPT
  if (piston != 0.0) 
    p_vel = 0.0;
#endif

  scale = 0.5 * time_step * time_step;
  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces_propagate_vel:\n",this_node));

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      /* Rescale forces: f_rescaled = 0.5*dt*dt * f_calculated */
      p[i].f.f[0] *= scale;
      p[i].f.f[1] *= scale;
      p[i].f.f[2] *= scale;

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: SCAL f = (%.3e,%.3e,%.3e) v_old = (%.3e,%.3e,%.3e)\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));

      for(j = 0; j < 3 ; j++) {
#ifdef EXTERNAL_FORCES
	if (!(p[i].l.ext_flag & COORD_FIXED(j)))
#endif
	  /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * f(t+dt) */
	  p[i].m.v[j] += p[i].f.f[j];
#ifdef NPT
	if (piston != 0.0) 
	  p_vel += SQR(p[i].m.v[j]);
#endif
#ifdef NEMD
	if(j==0) nemd_store_velocity(p[i]);
#endif
      }

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
    }
  }
#ifdef NPT
  finalize_p_inst_npt();
#endif

#ifdef NEMD
    nemd_change_momentum();
    nemd_store_velocity_profile();
#endif     

}

void finalize_p_inst_npt()
{
#ifdef NPT
  if (piston != 0.0) {
    double p_tmp=0.0;
    /* finalize derivation of p_inst */
    p_vel /= SQR(time_step);
    p_inst = p_vir+p_vel;
    MPI_Reduce(&p_inst, &p_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (this_node == 0) {
      p_inst = p_tmp/(3.0*NpT_volume);
      p_diff = p_diff  +  (p_inst-p_ext)*0.5*time_step + friction_thermo_nptiso();
// fprintf(stderr, "%d: %f -> %f \n",this_node,p_inst,p_diff);
    }
  }

#endif
}

void propagate_press_box_pos_and_rescale_npt()
{
#ifdef NPT
  if (piston != 0.0) {
    Cell *cell;
    Particle *p;
    int i, j, np, c;
    double scal[3], L_new;

    /* finalize derivation of p_inst */
    finalize_p_inst_npt();

    /* adjust NpT_volume; prepare pos- and vel-rescaling */
    if (this_node == 0) {
      NpT_volume += inv_piston*p_diff*0.5*time_step;
      scal[2] = SQR(box_l[0])/pow(NpT_volume,2.0/3.0);
      NpT_volume += inv_piston*p_diff*0.5*time_step;
      L_new = pow(NpT_volume,1.0/3.0);
// fprintf(stderr, "%d: %f -> %f; volume: t=%f => L=%f || ",this_node,p_inst,p_diff,NpT_volume,L_new);
      scal[1] = L_new/box_l[0];
      scal[0] = 1/scal[1];
      box_l[0] = box_l[1] = box_l[2] = L_new;
    }
    MPI_Bcast(scal,  3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    on_NpT_boxl_change();

    /* propagate positions while rescaling velocities */
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c]; p  = cell->part; np = cell->n;
      for(i = 0; i < np; i++) {	for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	  if (!(p[i].l.ext_flag & COORD_FIXED(j)))	
#endif
	    {
		p[i].r.p[j]  = scal[1]*(p[i].r.p[j] + scal[2]*p[i].m.v[j]);
		p[i].m.v[j] *= scal[0];
	    }
	}
	ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT:PV_1 v_new=(%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
	ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT:PPOS p=(%.3f,%.3f,%.3f)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2])); 
#ifdef ADDITIONAL_CHECKS
	force_and_velocity_check(&p[i]); 
#endif
	/* Verlet criterion check */
	if(distance2(p[i].r.p,p[i].l.p_old) > skin2 )
	  rebuild_verletlist = 1; 
      }
    }
  }
#endif
}

void propagate_vel_pos() 
{
  Cell *cell;
  Particle *p;
  int c, i, j, np;
#ifdef NPT
  p_vel = 0.0;
#endif

  INTEG_TRACE(fprintf(stderr,"%d: propagate_vel_pos:\n",this_node));
  rebuild_verletlist = 0;
  skin2 = SQR(skin/2.0);

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
#ifdef NEMD 
      nemd_add_velocity(&p[i]);
#endif
      for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	if (!(p[i].l.ext_flag & COORD_FIXED(j)))	
#endif
	  {
	    /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
	    p[i].m.v[j] += p[i].f.f[j];
#ifdef NPT
	    if (piston != 0.0) 
	      p_vel += SQR(p[i].m.v[j]);
	    else
#endif
	      /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
	      p[i].r.p[j] += p[i].m.v[j];
	  }
      }
 
      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PPOS p = (%.3f,%.3f,%.3f)\n",this_node,p[i].r.p[0],p[i].r.p[1],p[i].r.p[2])); 
      
#ifdef ADDITIONAL_CHECKS
      force_and_velocity_check(&p[i]);
#endif

      /* Verlet criterion check */
#ifdef NPT
      if(!(piston != 0.0))
#endif
	if(distance2(p[i].r.p,p[i].l.p_old) > skin2 )
	  rebuild_verletlist = 1; 
    }
  }

#ifdef NPT
  /* Propagate pressure, box_length (2 times) and positions, rescale
     positions and velocities and check verlet list criterion (only NPT) */
  propagate_press_box_pos_and_rescale_npt();
#endif


#ifdef ADDITIONAL_CHECKS
  force_and_velocity_display(); 
#endif

  /* communicate verlet criterion */
  anounce_rebuild_vlist();

#ifdef NEMD
  nemd_change_momentum() ;
#endif

}

void force_and_velocity_check(Particle *p) 
{
#ifdef ADDITIONAL_CHECKS
  double db_force,db_vel;
  /* force check */
  db_force = SQR(p->f.f[0])+SQR(p->f.f[1])+SQR(p->f.f[2]);
  if(db_force > skin2) 
    fprintf(stderr,"%d: Part %d has force %f (%f,%f,%f)\n",
	    this_node,p->p.identity,sqrt(db_force),
	    p->f.f[0],p->f.f[1],p->f.f[2]);
  if(db_force > db_max_force) { db_max_force=db_force; db_maxf_id=p->p.identity; }
  /* velocity check */
  db_vel   = SQR(p->m.v[0])+SQR(p->m.v[1])+SQR(p->m.v[2]);
  if(db_vel > skin2) 
	fprintf(stderr,"%d: Part %d has velocity %f (%f,%f,%f)\n",
		this_node,p->p.identity,sqrt(db_vel),
		p->m.v[0],p->m.v[1],p->m.v[2]);
  if(db_vel > db_max_vel) { db_max_vel=db_vel; db_maxv_id=p->p.identity; } 
#endif
}

void force_and_velocity_display() 
{
#ifdef ADDITIONAL_CHECKS
  if(db_max_force > skin2) 
    fprintf(stderr,"%d: max_force=%e, part=%d f=(%e,%e,%e)\n",this_node,
	    sqrt(db_max_force),db_maxf_id,local_particles[db_maxf_id]->f.f[0],
	    local_particles[db_maxf_id]->f.f[1],local_particles[db_maxf_id]->f.f[2]);
  if(db_max_vel > skin2)
    fprintf(stderr,"%d: max_vel=%e, part=%d v=(%e,%e,%e)\n",this_node,
	    sqrt(db_max_vel),db_maxv_id,local_particles[db_maxv_id]->m.v[0],
	    local_particles[db_maxv_id]->m.v[1],local_particles[db_maxv_id]->m.v[2]);
#endif
}
