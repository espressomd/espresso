/*****************************************************/
/*******************  INTEGRATE.C  *******************/
/*****************************************************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "integrate.h"
#include "communication.h"
#include "cells.h"
#include "verlet.h"
#include "ghosts.h"
#include "forces.h"
#include "debug.h"
#include "p3m.h"
#include "utils.h"

/*******************  variables  *******************/

double time_step = 0.001;
double max_cut = 2.0;
double skin = 0.4;
double max_range;
double max_range2;
int calc_forces_first = 1;
/* extern p3m_struct p3m; */

/*******************  functions  *******************/
void rescale_forces();
void propagate_velocities();
void propagate_positions(); 

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
  if (!strncmp(argv[1], "init", strlen(argv[1]))) n_steps=0;
  else if (!strncmp(argv[1], "exit", strlen(argv[1]))) n_steps=-1;
  else n_steps = atol(argv[1]);

  if(n_steps < -1) {
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
  int i;

  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_init\n",this_node));
  INTEG_TRACE(fprintf(stderr,"%d: n_node =%d npart=%d\n",
		      this_node,n_nodes,max_seen_particle));
  INTEG_TRACE(fprintf(stderr,"%d: n_particles = %d\n",
		      this_node, n_particles));
  max_range  = max_cut + skin;
  max_range2 = max_range* max_range;

  /* initialize link cell structure */
  cells_init();  fflush(stderr); MPI_Barrier(MPI_COMM_WORLD);
  /* allocate and initialize local indizes */
  local_index = (int *)malloc((max_seen_particle + 1)*sizeof(int));
  for(i=0;i<=max_seen_particle;i++) local_index[i] = -1;
  for(i=0;i<n_particles;i++) local_index[particles[i].identity] = i;
  /* initialize ghost structure */
  ghost_init(); fflush(stderr);  MPI_Barrier(MPI_COMM_WORLD);
  /* initialize verlet list structure */
  verlet_init(); fflush(stderr);  MPI_Barrier(MPI_COMM_WORLD);
  /* initialize force structure */
  force_init(); fflush(stderr);  MPI_Barrier(MPI_COMM_WORLD);
  /* initialize p3m */
  P3M_init();
}

void integrate_vv(int n_steps)
{
  int i;
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv: %d steps\n",this_node,
		      n_steps));
  /* this is just for security reasons */
  MPI_Barrier(MPI_COMM_WORLD);
  /* check init */
  if(rebuild_verletlist == 1) {
    exchange_part();
    sort_particles_into_cells();
    exchange_ghost();
    build_verlet_list();
  }
  if(calc_forces_first == 1) {
    force_calc();
    collect_ghost_forces();
    rescale_forces();
    calc_forces_first = 0;
  }

  /* integration loop */
  INTEG_TRACE(printf("%d START INTEGRATION\n",this_node));
  for(i=0;i<n_steps;i++) {
    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,i));
    propagate_velocities();
    propagate_positions();
   /* rebuild_verletlist = 1; */
    if(rebuild_verletlist == 1) {
      exchange_part();
      sort_particles_into_cells();
      exchange_ghost();
      build_verlet_list();
    }
    else {
      update_ghost_pos();
    }
    force_calc();
    collect_ghost_forces();
    rescale_forces();
    propagate_velocities();
  }

}

void integrate_vv_exit()
{
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_exit\n",this_node));
  cells_exit();
  free(local_index);
  ghost_exit();
  verlet_exit();
  force_exit();
}

void rescale_forces()
{
  int i,d;
  double scale;
  scale = 0.5 * time_step * time_step;
  INTEG_TRACE(fprintf(stderr,"%d: rescale_forces:\n",this_node));
  for(i=0;i<n_particles;i++)
    for(d=0; d<3; d++) particles[i].f[d] *= scale;
}

void propagate_velocities() 
{
  int i;
  INTEG_TRACE(fprintf(stderr,"%d: propagate_velocities:\n",this_node));
  for(i=0;i<n_particles;i++)
    {  
      particles[i].v[0] += particles[i].f[0];
      particles[i].v[1] += particles[i].f[1];
      particles[i].v[2] += particles[i].f[2];
    }
}

void propagate_positions() 
{
  int i;
  int *verlet_flags = malloc(sizeof(int)*n_nodes);
  double d2[3], dist2, skin2;
  INTEG_TRACE(fprintf(stderr,"%d: propagate_positions:\n",this_node));
  rebuild_verletlist = 0;
  skin2 = SQR(skin);

  for(i=0;i<n_particles;i++)
    {  
      particles[i].p[0] += particles[i].v[0];
      particles[i].p[1] += particles[i].v[1];
      particles[i].p[2] += particles[i].v[2];
      d2[0] = SQR(particles[i].p[0] - particles[i].p_old[0]);
      d2[1] = SQR(particles[i].p[1] - particles[i].p_old[1]);
      d2[2] = SQR(particles[i].p[2] - particles[i].p_old[2]);
      dist2 = d2[0] + d2[1] + d2[2];
      INTEG_TRACE(fprintf(stderr,"%d: prop_pos: P[%d] moved %f\n",this_node,i,sqrt(dist2)));
      if( dist2 > skin2 )
	rebuild_verletlist = 1; 
    }

  MPI_Gather(&rebuild_verletlist, 1, MPI_INT, verlet_flags, 1, MPI_INT, 0, MPI_COMM_WORLD);
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

