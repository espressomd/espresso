// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
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

int    integ_switch     = INTEG_METHOD_NVT;

int n_verlet_updates    = 0;

double time_step        = -1.0;
double sim_time         = 0.0;
double skin             = -1.0;
double skin2;
double max_range        = -1.0;
double max_range2       = -1.0;

int    resort_particles = 1;
int    recalc_forces    = 1;

double verlet_reuse     = 0.0;

#ifdef ADDITIONAL_CHECKS
double db_max_force = 0.0, db_max_vel = 0.0;
int    db_maxf_id   = 0,   db_maxv_id = 0;
#endif

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Rescale all particle forces with \f[ 0.5 \Delta t^2 \f]. */
void rescale_forces();
/** Propagate the velocities. Integration step 1 of the Velocity Verlet integrator:<br>
    \f[ v(t+0.5 \Delta t) = v(t) + 0.5 \Delta t f(t) \f] */
void propagate_vel();
/** Propagate the positions. Integration step 2 of the Velocity Verletintegrator:<br>
    \f[ p(t+\Delta t) = p(t) + \Delta t  v(t+0.5 \Delta t) \f] */
void propagate_pos();
/** Propagate the velocities and positions. Integration step 1 and 2
    of the Velocity Verlet integrator: <br>
    \f[ v(t+0.5 \Delta t) = v(t) + 0.5 \Delta t f(t) \f] <br>
    \f[ p(t+\Delta t) = p(t) + \Delta t  v(t+0.5 \Delta t) \f] */
void propagate_vel_pos();
/** Rescale all particle forces with \f[ 0.5 \Delta t^2 \f] and propagate the velocities.
    Integration step 4 of the Velocity Verletintegrator:<br> 
    \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t) \f] */
void rescale_forces_propagate_vel();

/** Integrator stability check (see compile flag ADDITIONAL_CHECKS). */
void force_and_velocity_check(Particle *p); 
/** Integrator stability check (see compile flag ADDITIONAL_CHECKS). */
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

/**  Hand over integrate usage information to tcl interpreter. */
int integrate_usage(Tcl_Interp *interp) 
{
  Tcl_AppendResult(interp, "Usage of tcl-command integrate:\n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate <INT n steps>' for integrating n steps \n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate set' for printing integrator status \n", (char *)NULL);
  Tcl_AppendResult(interp, "'integrate set nvt' for enabling NVT integration or \n" , (char *)NULL);
#ifdef NPT
  Tcl_AppendResult(interp, "'integrate set npt_isotropic <DOUBLE p_ext> [<DOUBLE piston>]' for enabling isotropic NPT integration \n" , (char *)NULL);
#endif
  return (TCL_OK);
}

/** Hand over integrate status information to tcl interpreter. */
int integrate_print_status(Tcl_Interp *interp) 
{
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE];
  switch (integ_switch) {
  case INTEG_METHOD_NVT:
    Tcl_AppendResult(interp, "{ set nvt }", (char *)NULL);
    return (TCL_OK);
  case INTEG_METHOD_NPT_ISO:
    Tcl_PrintDouble(interp, nptiso.p_ext, buffer);
    Tcl_AppendResult(interp, "{ set npt_isotropic ", buffer, (char *)NULL);
    Tcl_PrintDouble(interp, nptiso.piston, buffer);
    Tcl_AppendResult(interp, " ", buffer, " } ", (char *)NULL);
    return (TCL_OK);
  }
  return (TCL_ERROR);
}

/** Parse integrate nvt command */
int integrate_parse_nvt(Tcl_Interp *interp, int argc, char **argv)
{
  integ_switch = INTEG_METHOD_NVT;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
  return (TCL_OK);
}

/** Parse integrate npt_isotropic command */
int integrate_parse_npt_isotropic(Tcl_Interp *interp, int argc, char **argv)
{
  if (argc < 4) {
    Tcl_AppendResult(interp, "wrong # args: \n", (char *)NULL);
    return integrate_usage(interp);
  }  
  /* set parameters p_ext and piston */
  if ( !ARG_IS_D(3, nptiso.p_ext) )  return integrate_usage(interp);
  if ( argc > 4 ) { 
    if(!ARG_IS_D(4, nptiso.piston) ) return integrate_usage(interp);
    piston_callback(interp, &nptiso.piston); }
  else if ( nptiso.piston <= 0.0 ) {
    Tcl_AppendResult(interp, "You must set <piston> as well before you can use this integrator! \n", (char *)NULL);
    return integrate_usage(interp);
  }
  p_ext_callback(interp, &nptiso.p_ext);
  /* set integrator switch */
  integ_switch = INTEG_METHOD_NPT_ISO;
  mpi_bcast_parameter(FIELD_INTEG_SWITCH);
  return (TCL_OK);
}

int integrate(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  int  n_steps;
  
  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if (argc < 1) {
    Tcl_AppendResult(interp, "wrong # args: \n\"", (char *) NULL);
    return integrate_usage(interp);  }
  else if (argc < 2) {                    return integrate_print_status(interp); }

  if (ARG1_IS_S("set")) {
    if      (argc < 3)                    return integrate_print_status(interp);
    if      (ARG_IS_S(2,"nvt"))           return integrate_parse_nvt(interp, argc, argv);
#ifdef NPT
    else if (ARG_IS_S(2,"npt_isotropic")) return integrate_parse_npt_isotropic(interp, argc, argv);
#endif
    else {
      Tcl_AppendResult(interp, "unknown integrator method:\n", (char *)NULL);
      return integrate_usage(interp);
    }
  } else if ( !ARG_IS_I(1,n_steps) ) return integrate_usage(interp);

  /* go on with integrate <n_steps> */
  if(n_steps < 0) {
    Tcl_AppendResult(interp, "illegal number of steps\n", (char *) NULL);
    return integrate_usage(interp);;
  }
  /* perform integration */
  mpi_integrate(n_steps);

  return (TCL_OK);
}

/************************************************************/

void integrate_vv_recalc_maxrange()
{
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_recalc_maxrange:\n",this_node));

  /* maximal interaction cutoff */
  calc_maximal_cutoff();
  if (max_cut < 0.0) {
    max_range  = -1.0;
    max_range2 = -1.0;
    return;
  }
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
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    /* prepare NpT-integration */
    nptiso.inv_piston = 1/(1.0*nptiso.piston);
    nptiso.volume     = box_l[0]*box_l[1]*box_l[2];
    if (recalc_forces) { 
      nptiso.p_inst = 0.0;  
      nptiso.p_vir  = 0.0; 
      nptiso.p_vel  = 0.0; 
    }
  }
#endif
}

/************************************************************/

void integrate_vv(int n_steps)
{
  int i;

  /* Prepare the Integrator */
  on_integration_start();
  /* Verlet list criterion */
  skin2 = SQR(skin/2.0);

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

  n_verlet_updates = 0;

  /* Integration loop */
  for(i=0;i<n_steps;i++) {

    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,i));

    /* Integration Steps: Step 1 and 2 of Velocity Verlet scheme:
       v(t+0.5*dt) = v(t) + 0.5*dt * f(t)
       p(t + dt)   = p(t) + dt * v(t+0.5*dt)
       NOTE 1: Prefactors do not occur in formulas since we use 
               rescaled forces and velocities. 
       NOTE 2: Depending on the integration method Step 1 and Step 2 
               can not be combined for the translation. 
    */
    if(integ_switch == INTEG_METHOD_NPT_ISO || nemd_method != NEMD_METHOD_OFF) {
      propagate_vel();  propagate_pos(); }
    else
      propagate_vel_pos();
#ifdef ROTATION
    propagate_omega_quat(); 
#endif

    cells_update_ghosts();

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
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    MPI_Bcast(&nptiso.p_inst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nptiso.p_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nptiso.volume, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
  if(integ_switch == INTEG_METHOD_NPT_ISO)
    nptiso.p_vel = 0.0;
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
#ifdef NPT
	  if(integ_switch == INTEG_METHOD_NPT_ISO) {
	    nptiso.p_vel += SQR(p[i].m.v[j]);
	    p[i].m.v[j] += p[i].f.f[j] + friction_therm0_nptiso(p[i].m.v[j]); }
	  else
#endif
	    /* Propagate velocity: v(t+dt) = v(t+0.5*dt) + 0.5*dt * f(t+dt) */
	    p[i].m.v[j] += p[i].f.f[j];
      }

      ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_2 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
    }
  }
#ifdef NPT
  finalize_p_inst_npt();
  // if((integ_switch == INTEG_METHOD_NPT_ISO) && (this_node==0)) fprintf(stderr,"%d/B: p_inst=%f \n",this_node,nptiso.p_inst);
#endif
}

void finalize_p_inst_npt()
{
#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    double p_tmp=0.0;
    /* finalize derivation of p_inst */
    nptiso.p_vel /= SQR(time_step);
    nptiso.p_inst = nptiso.p_vir + nptiso.p_vel;
    // fprintf(stderr,"%d: p_inst=%f+%f=%f",this_node,nptiso.p_vir/(3.*nptiso.volume),nptiso.p_vel/(3.*nptiso.volume),nptiso.p_inst/(3.*nptiso.volume));
    MPI_Reduce(&nptiso.p_inst, &p_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (this_node == 0) {
      nptiso.p_inst = p_tmp/(3.0*nptiso.volume);
      nptiso.p_diff = nptiso.p_diff  +  (nptiso.p_inst-nptiso.p_ext)*0.5*time_step + friction_thermV_nptiso(nptiso.p_diff);
      // fprintf(stderr, " = %f -> %f  (volume: %f)\n",nptiso.p_inst,nptiso.p_diff,nptiso.volume/(box_l[0]*box_l[1]*box_l[2]));
    }
  }

#endif
}

void propagate_press_box_pos_and_rescale_npt()
{
#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    Cell *cell;
    Particle *p;
    int i, j, np, c;
    double scal[3], L_new;

    /* finalize derivation of p_inst */
    // fprintf(stderr, "%d: propagate_press_box...: finalize_p_inst ",this_node);
    finalize_p_inst_npt();

    /* adjust \ref nptiso_struct::nptiso.volume; prepare pos- and vel-rescaling */
    if (this_node == 0) {
      nptiso.volume += nptiso.inv_piston*nptiso.p_diff*0.5*time_step;
      scal[2] = SQR(box_l[0])/pow(nptiso.volume,2.0/3.0);
      nptiso.volume += nptiso.inv_piston*nptiso.p_diff*0.5*time_step;
      L_new = pow(nptiso.volume,1.0/3.0);
      if (nptiso.volume < 0.0) { 
	fprintf(stderr, "%d: ERROR: Your choice of piston=%f, dt=%f, p_diff=%f just caused the volume to become negative!\nTry decreasing dt...\n",\
		this_node,nptiso.piston,time_step,nptiso.p_diff); errexit(); }
      //else fprintf(stderr, "Q=%f, dt=%f, p_diff=%f => L=%f (V=%f; p_inst=%f)\n",nptiso.piston,time_step,nptiso.p_diff,L_new,nptiso.volume,nptiso.p_inst);
      scal[1] = L_new/box_l[0];
      scal[0] = 1/scal[1];
      //box_l[0] = box_l[1] = box_l[2] = L_new;
    }
    //    scal[0] = scal[1] = scal[2] = 1.0;
    //    nptiso.volume = box_l[0]*box_l[1]*box_l[2];
    MPI_Bcast(scal,  3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //on_NpT_boxl_change(scal[1]);

    /* propagate positions while rescaling positions and velocities */
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c]; p  = cell->part; np = cell->n;
      for(i = 0; i < np; i++) {	for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	  if (!(p[i].l.ext_flag & COORD_FIXED(j)))	
#endif
	    {
		p[i].r.p[j]      = scal[1]*(p[i].r.p[j] + scal[2]*p[i].m.v[j]);
		p[i].l.p_old[i] *= scal[1];
		p[i].m.v[j]     *= scal[0];
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

    /* Apply new volume to the box-length, communicate it, and account for necessary adjustments to the cell geometry */
    if (this_node == 0) box_l[0] = box_l[1] = box_l[2] = L_new;
    MPI_Bcast(box_l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //    on_NpT_boxl_change(scal[1]);
    on_NpT_boxl_change(0.0);
  }
#endif
}

void propagate_vel()
{
  Cell *cell;
  Particle *p;
  int c, i, j, np;
#ifdef NPT
  nptiso.p_vel = 0.0;
#endif

  INTEG_TRACE(fprintf(stderr,"%d: propagate_vel:\n",this_node));

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	if (!(p[i].l.ext_flag & COORD_FIXED(j)))	
#endif
	  {
#ifdef NPT
	    if(integ_switch == INTEG_METHOD_NPT_ISO) {
	      p[i].m.v[j] += p[i].f.f[j] + friction_therm0_nptiso(p[i].m.v[j]);
	      nptiso.p_vel += SQR(p[i].m.v[j]); }
	    else
#endif
	      /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
	      p[i].m.v[j] += p[i].f.f[j];

	    /* SPECIAL TASKS in particle loop */
#ifdef NEMD
	    if(j==0) nemd_get_velocity(p[i]);
#endif
	  }

	ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: PV_1 v_new = (%.3e,%.3e,%.3e)\n",this_node,p[i].m.v[0],p[i].m.v[1],p[i].m.v[2]));
#ifdef ADDITIONAL_CHECKS
      force_and_velocity_check(&p[i]);
#endif
      }
    }
  }
#ifdef ADDITIONAL_CHECKS
  force_and_velocity_display(); 
#endif

  /* SPECIAL TASKS after velocity propagation */
#ifdef NEMD
  nemd_change_momentum();
  nemd_store_velocity_profile();
#endif
}

void propagate_pos() 
{
  INTEG_TRACE(fprintf(stderr,"%d: propagate_pos:\n",this_node));
  if(integ_switch == INTEG_METHOD_NPT_ISO) 
    /* Special propagator for NPT ISOTROPIC */
    /* Propagate pressure, box_length (2 times) and positions, rescale
       positions and velocities and check verlet list criterion (only NPT) */
    propagate_press_box_pos_and_rescale_npt();
  else {
    Cell *cell;
    Particle *p;
    int c, i, j, np;
    
    rebuild_verletlist = 0;

    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i = 0; i < np; i++) {
	for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	  if (!(p[i].l.ext_flag & COORD_FIXED(j)))	
#endif
	    {
#ifdef NEMD
	      /* change momentum of each particle in top and bottom slab */
	      if(j==0) nemd_add_velocity(&p[i]);
#endif	    
	      /* Propagate positions (only NVT): p(t + dt)   = p(t) + dt * v(t+0.5*dt) */
	      p[i].r.p[j] += p[i].m.v[j];
	    }
	}
	/* Verlet criterion check */
	if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) rebuild_verletlist = 1; 
      }
    }
  }
  /* communicate verlet criterion */
  anounce_rebuild_vlist();
}

void propagate_vel_pos() 
{
  Cell *cell;
  Particle *p;
  int c, i, j, np;

  INTEG_TRACE(fprintf(stderr,"%d: propagate_vel_pos:\n",this_node));
  rebuild_verletlist = 0;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      for(j=0; j < 3; j++){
#ifdef EXTERNAL_FORCES
	if (!(p[i].l.ext_flag & COORD_FIXED(j)))	
#endif
	  {
	    /* Propagate velocities: v(t+0.5*dt) = v(t) + 0.5*dt * f(t) */
	    p[i].m.v[j] += p[i].f.f[j];

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
      if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) rebuild_verletlist = 1; 
    }
  }

#ifdef ADDITIONAL_CHECKS
  force_and_velocity_display();
#endif

  /* communicate verlet criterion */
  anounce_rebuild_vlist();
}

void force_and_velocity_check(Particle *p) 
{
#ifdef ADDITIONAL_CHECKS
  int i;
  double db_force,db_vel;
  /* distance_check */
  for (i = 0; i < 3; i++)
    if(fabs(p->r.p[i] - p->l.p_old[i]) > local_box_l[i]) {
      /* put any cellsystems here which do not rely on rebuild_verletlist */
      if (cell_structure.type != CELL_STRUCTURE_NSQUARE) {
	fprintf(stderr, "%d: particle %d moved further than local box length by %lf %lf %lf\n",
		this_node, p->p.identity, p->r.p[0] - p->l.p_old[0], p->r.p[1] - p->l.p_old[1],
		p->r.p[2] - p->l.p_old[2]);
      }
    }

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
