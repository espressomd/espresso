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
#include "communication.h"
#include "grid.h"
#include "cells.h"
#include "verlet.h"
#include "ghosts.h"
#include "forces.h"
#include "debug.h"
#include "p3m.h"
#include "utils.h"
#include "thermostat.h"

/************************************************
 * DEFINES
 ************************************************/

/* MPI tags for the fft communications: */
/** Tag for communication in verlet fix: propagate_positions()  */
#define REQ_INT_VERLET   400

/*******************  variables  *******************/

double time_step = -1.0;
double skin = -1.0;
double max_range;
double max_range2;
int    calc_forces_first;
int    vv_integrator_initialized = 0;

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Rescale all particle forces with \f[ 0.5 \Delta t^2 \f]. */
void rescale_forces();
/** Propagate velocities: \f[ {\bf v} += {\bf f} \f]. */
void propagate_velocities();
/** Propagate positions: \f[ {\bf pos} += {\bf v} \f]. 
    Here also the verlet criterion is checked and communicated.
    \todo Put verlet criterion as inline function to verlet.h if possible.
*/
void propagate_positions(); 

void propagate_vel_pos();
void rescale_forces_propagate_vel(); 


void print_local_index();

/*@}*/
/************************************************************/

int integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv)
{
  int  n_steps;
  
  INTEG_TRACE(fprintf(stderr,"%d: integrate:\n",this_node));

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <step num> \"", (char *) NULL);
    return (TCL_ERROR);
  }

  /* translate argument */
  if (!strncmp(argv[1], "init", strlen(argv[1]))) n_steps=-1;
  else if (!strncmp(argv[1], "exit", strlen(argv[1]))) n_steps=-2;
  else n_steps = atol(argv[1]);

  if(n_steps < -2) {
    Tcl_AppendResult(interp, "illegal number of steps", (char *) NULL);
    return (TCL_ERROR);
  }

  /* flush remaining information from master node to slaves */
  particle_finalize_data();

  /* assume velocity verlet integration with langevin thermostat */
  if (argc != 2) {
    Tcl_AppendResult(interp, "too many arguments:  should be \"",
		     argv[0], " <task> \"", (char *) NULL);
    return (TCL_ERROR);
  }

  mpi_integrate(n_steps);

  return (TCL_OK);
}

void integrate_vv_init()
{
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_init:\n",this_node));
  
  /* sanity checks */
  if(time_step < 0.0 || skin < 0.0 ) {
    fprintf(stderr,"%d: ERROR: Can not initialize the integrator!:\n",this_node);
    if( time_step < 0.0 )
      fprintf(stderr,"%d: PROBLEM: You have to set the time_step!\n",this_node);
    if( skin < 0.0 )
      fprintf(stderr,"%d: PROBLEM: You have to set the skin!\n",this_node);
    errexit();
  }
  else {
    
    /* maximal interaction cutoff */
    calc_maximal_cutoff();
    calc_minimal_box_dimensions();
    max_range  = max_cut + skin;
    max_range2 = max_range* max_range;

   /* check real space interaction cutoff/range*/
    if(max_cut < 0.0) {
      fprintf(stderr,"%d: ERROR: You have to specify at least one interaction\n",this_node);
      errexit();
    }
    if((min_box_l/2.0) < max_range) {
      fprintf(stderr,"%d: ERROR: Maximal real space interaction %f is larger than half of the minimal box dimension %f\n",this_node,max_range,min_box_l);
      errexit();
    }

    if(min_local_box_l < max_range) {
      fprintf(stderr,"%d: ERROR: Maximal real space interaction %f is larger than minimal local box length %f\n",this_node,max_range,min_local_box_l);
      errexit();
    }

#ifdef INTEG_DEBUG 
    if(this_node==0) {
      double n_inter;
      n_inter = (4.0/3.0)* PI * (max_range*max_range*max_range); /* interaction volume */
      n_inter /= (box_l[0]*box_l[1]*box_l[2]); /* particle density */
      n_inter /= 2.0;
      fprintf(stderr,"0: max_cut=%f, max_range+%f, n_inter/np^2=%f\n",max_cut,max_range,n_inter);
      fflush(stderr);
    }
#endif

    /* start initialization */

    /* initialize link cell structure */
    cells_re_init();  
    /* allocate and initialize local indizes */
    local_particles_init();
    /* initialize ghost structure */
    ghost_init(); 
    /* initialize force structure */
    force_init(); 
    /* initialize p3m */
    P3M_init();
    thermo_init();
    /* update integrator status */
    calc_forces_first         = 1;
    rebuild_verletlist        = 1;
    vv_integrator_initialized = 1;
  }
}

void integrate_vv(int n_steps)
{
  int i;

  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv:\n",this_node));
  if(vv_integrator_initialized == 0) {
    fprintf(stderr,"%d: ERROR: Velocity verlet Integrator not initialized\n",
	    this_node);
    errexit();
  }
  else {
    INTEG_TRACE(fprintf(stderr,"%d: integrate_vv: integrating %d steps\n",this_node,
			n_steps));
    /* check init */
    if(rebuild_verletlist == 1) {
      /* already here since exchange_part may destroy the ghost particles */
      invalidate_ghosts();
      exchange_part();
      sort_particles_into_cells();
      exchange_ghost();
      build_verlet_lists();
    }
    if(calc_forces_first == 1) {
      force_calc();
      collect_ghost_forces();
      rescale_forces();
      calc_forces_first = 0;
    }
    
    /* integration loop */
    INTEG_TRACE(fprintf(stderr,"%d START INTEGRATION\n",this_node));
    for(i=0;i<n_steps;i++) {
      INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,i));
      propagate_vel_pos();
      if(rebuild_verletlist == 1) {
	INTEG_TRACE(fprintf(stderr,"%d: Rebuild Verlet List\n",this_node));
	invalidate_ghosts();
	exchange_part();
	sort_particles_into_cells(); 
	exchange_ghost();
	build_verlet_lists();
      }
      else {
	update_ghost_pos();
      }
      force_calc();
      collect_ghost_forces();
      rescale_forces_propagate_vel();
    }
  }
}

void integrate_vv_exit()
{
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_exit\n",this_node));
  local_particles_exit();
  cells_exit();
  ghost_exit();
  force_exit();
  vv_integrator_initialized = 0;
}

/* Callback functions */
/************************************************************/


int skin_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "skin must be positiv.", (char *) NULL);
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
    Tcl_AppendResult(interp, "time step must be positiv.", (char *) NULL);
    return (TCL_ERROR);
  }
  time_step = data;
  mpi_bcast_parameter(FIELD_TIME_STEP);
  return (TCL_OK);
}

int calc_forces_first_callback(Tcl_Interp *interp, void *_data)
{
  int data = *(int *)_data;
  if (data != 0 && data != 1) {
    Tcl_AppendResult(interp, "calc_forces_first is an integrator flag and must be 0 or 1.", (char *) NULL);
    return (TCL_ERROR);
  }
  calc_forces_first = data;
  mpi_bcast_parameter(FIELD_CALC_FORCES_FIRST);
  return (TCL_OK);
}

/************************************************************/

void rescale_forces()
{
  Particle *p;
  int m,n,o,i, np;
  double scale;
  scale = 0.5 * time_step * time_step;
  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces:\n",this_node));
  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
      p[i].f[0] *= scale;
      p[i].f[1] *= scale;
      p[i].f[2] *= scale;
    }
  }
}

void propagate_velocities() 
{
  Particle *p;
  int m,n,o,i, np;
  INTEG_TRACE(fprintf(stderr,"%d: propagate_velocities:\n",this_node));
  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
      p[i].v[0] += p[i].f[0];
      p[i].v[1] += p[i].f[1];
      p[i].v[2] += p[i].f[2];
    }
  }
}

void rescale_forces_propagate_vel() 
{
  Particle *p;
  int m,n,o,i, np;
  double scale;

  scale = 0.5 * time_step * time_step;
  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces_propagate_vel:\n",this_node));

  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
      p[i].f[0] *= scale;
      p[i].f[1] *= scale;
      p[i].f[2] *= scale;

      p[i].v[0] += p[i].f[0];
      p[i].v[1] += p[i].f[1];
      p[i].v[2] += p[i].f[2];
    }
  }
}

void propagate_positions() 
{
  Particle *p;
  int m,n,o,i, np;

  int *verlet_flags = malloc(sizeof(int)*n_nodes);
  double skin2;
  INTEG_TRACE(fprintf(stderr,"%d: propagate_positions:\n",this_node));
  rebuild_verletlist = 0;
  skin2 = SQR(skin);

  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {
      p[i].r.p[0] += p[i].v[0];
      p[i].r.p[1] += p[i].v[1];
      p[i].r.p[2] += p[i].v[2];

      /* Verlet criterion check */
      if(distance2(p[i].r.p,p[i].p_old) > skin2 )
	rebuild_verletlist = 1; 
    }
  }

  /* communicate verlet criterion */
  MPI_Gather(&rebuild_verletlist, 1, MPI_INT, verlet_flags, 1, 
	     MPI_INT, 0, MPI_COMM_WORLD);
  if(this_node == 0)
    {
      rebuild_verletlist = 0;
      for(i=0;i<n_nodes;i++)
	if(verlet_flags[i]>0)
	  {
	    rebuild_verletlist = 1;
	    break;
	  }
    }
  MPI_Barrier(MPI_COMM_WORLD); 
  MPI_Bcast(&rebuild_verletlist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  INTEG_TRACE(fprintf(stderr,"%d: prop_pos: rebuild_verletlist=%d\n",this_node,rebuild_verletlist));

  free(verlet_flags);
}

void propagate_vel_pos() 
{
  Particle *p;
  int m,n,o,i, np;
  int *verlet_flags = malloc(sizeof(int)*n_nodes);
  double skin2;

#ifdef ADDITIONAL_CHECKS
  double db_force,db_vel,e_kin=0.0,tot_kin_energy=0.0;
  double db_max_force=0.0, db_max_vel=0.0;
  int db_maxf_id=0,db_maxv_id=0;
  double *kin_energies=NULL;
  FILE *stat;
#endif

  INTEG_TRACE(fprintf(stderr,"%d: propagate_vel_pos:\n",this_node));
  rebuild_verletlist = 0;
  skin2 = SQR(skin);

  INNER_CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    for(i = 0; i < np; i++) {

      p[i].v[0] += p[i].f[0];
      p[i].v[1] += p[i].f[1];
      p[i].v[2] += p[i].f[2];

#ifdef ADDITIONAL_CHECKS
      /* fprintf(stderr,"%d: Propagate P %d from (%.4f,%.4f,%.4f) to (%.4f,%.4f,%.4f)\n",this_node,
	 p[i].r.identity, p[i].r.p[0], p[i].r.p[1], p[i].r.p[2],
	 p[i].r.p[0]+p[i].v[0], p[i].r.p[1]+p[i].v[1], p[i].r.p[2]+p[i].v[2]); */
#endif

      p[i].r.p[0] += p[i].v[0];
      p[i].r.p[1] += p[i].v[1];
      p[i].r.p[2] += p[i].v[2];

#ifdef ADDITIONAL_CHECKS
      /* force check */
      db_force = SQR(p[i].f[0])+SQR(p[i].f[1])+SQR(p[i].f[2]);
      if(sqrt(db_force) > skin) 
	fprintf(stderr,"%d: Part %d has force %f (%f,%f,%f)\n",
		this_node,p[i].r.identity,sqrt(db_force),
		p[i].f[0],p[i].f[1],p[i].f[2]);
      if(db_force > db_max_force) { db_max_force=db_force; db_maxf_id=p[i].r.identity; }
      /* velocity check */
      db_vel   = SQR(p[i].v[0])+SQR(p[i].v[1])+SQR(p[i].v[2]);
      if(sqrt(db_vel) > skin) 
	fprintf(stderr,"%d: Part %d has velocity %f (%f,%f,%f)\n",
		this_node,p[i].r.identity,sqrt(db_vel),
		p[i].v[0],p[i].v[1],p[i].v[2]);
      if(db_vel > db_max_vel) { db_max_vel=db_vel; db_maxv_id=p[i].r.identity; }
      /* kinetic energy */
      e_kin += db_vel;
#endif

      /* Verlet criterion check */
      if(distance2(p[i].r.p,p[i].p_old) > skin2 )
	rebuild_verletlist = 1; 
    }
  }

#ifdef ADDITIONAL_CHECKS
  if(this_node==0) {
    stat=fopen("e_kin.dat","a");
    kin_energies = malloc(sizeof(double)*n_nodes);
    MPI_Gather(&e_kin, 1, MPI_DOUBLE, kin_energies, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(i=0;i<n_nodes;i++) tot_kin_energy += kin_energies[i];
    fprintf(stat,"%e\n",tot_kin_energy/(2.0*time_step*time_step));
    free(kin_energies);
    fclose(stat); 
  }
  else {
    MPI_Gather(&e_kin, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  if(sqrt(db_max_force)>1e-5) 
    fprintf(stderr,"%d: max_force=%e, part=%d f=(%e,%e,%e)\n",this_node,
	    sqrt(db_max_force),db_maxf_id,local_particles[db_maxf_id]->f[0],
	    local_particles[db_maxf_id]->f[1],local_particles[db_maxf_id]->f[2]);
  if(sqrt(db_max_vel)>1e-3)
    fprintf(stderr,"%d: max_vel=%e, part=%d v=(%e,%e,%e)\n",this_node,
	    sqrt(db_max_vel),db_maxv_id,local_particles[db_maxv_id]->v[0],
	    local_particles[db_maxv_id]->v[1],local_particles[db_maxv_id]->v[2]);
#endif

  /* communicate verlet criterion */
  MPI_Gather(&rebuild_verletlist, 1, MPI_INT, verlet_flags, 1, 
	     MPI_INT, 0, MPI_COMM_WORLD);
  if(this_node == 0)
    {
      rebuild_verletlist = 0;
      for(i=0;i<n_nodes;i++)
	if(verlet_flags[i]>0)
	  {
	    rebuild_verletlist = 1;
	    break;
	  }
    }
  MPI_Barrier(MPI_COMM_WORLD); 
  MPI_Bcast(&rebuild_verletlist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  INTEG_TRACE(fprintf(stderr,"%d: prop_pos: rebuild_verletlist=%d\n",this_node,rebuild_verletlist));

  free(verlet_flags);
}
