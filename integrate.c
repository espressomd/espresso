/*****************************************************/
/*******************  INTEGRATE.C  *******************/
/*****************************************************/

#include "integrate.h"

//#define DEBUG

/*******************  variables  *******************/

double time_step = 0.001;
double max_cut = 2.0;
double skin = 0.4;
double max_range;
double max_range2;

int calc_forces_first = 1;

/*******************  functions  *******************/

void propagate_velocities();
void propagate_positions(); 

void tcl_integrator_init(Tcl_Interp *interp)
{
  Tcl_CreateCommand(interp, "integrate", integrate, 0, NULL);
}

int integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv)
{
  int  n_steps;
  char buffer[256];
  
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: integrate:\n",this_node); 
#endif

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <step num> \"", (char *) NULL);
    return (TCL_ERROR);
  }

  n_steps = atol(argv[1]);
  if(n_steps < -1) {
    Tcl_AppendResult(interp, "illegal number of steps", (char *) NULL);
    return (TCL_ERROR);
  }

  /* assume velocity verlet integration with langevin thermostat */
  if (argc == 2) {
    mpi_integrate(n_steps);
    return (TCL_OK);
  }
  else {
    Tcl_AppendResult(interp, "too many arguments:  should be \"",
		     argv[0], " <step num> \"", (char *) NULL);
    return (TCL_ERROR);
  }
}


void integrate_vv_init()
{
  int i;

#ifdef DEBUG
  if(this_node < 2) fprintf(stderr,"%d: integrate_vv_init\n",this_node);
  if(this_node < 2) fprintf(stderr,"%d: nproc =%d npart=%d\n",
			   this_node,nprocs,n_total_particles);
  fprintf(stderr,"%d: n_total_part=%d  n_particles = %d\n",
			   this_node,n_total_particles,n_particles);
#endif
  max_range  = max_cut + skin;
  max_range2 = max_range* max_range;

  /* initialize link cell structure */
  cells_init();
  /* allocate and initialize local indizes */
  local_index = (int *)malloc(n_total_particles*sizeof(int));
  for(i=0;i<n_total_particles;i++) local_index[i] = -1;
  for(i=0;i<n_particles;i++) local_index[particles[i].identity] = i;
  /* initialize ghost structure */
  ghost_init();
  /* initialize verlet list structure */
  verlet_init();
  /* initialize force structure */
  force_init();
}

void integrate_vv(int n_steps)
{
  int i,j;
#ifdef DEBUG
  if(this_node<2) fprintf(stderr,"%d: integrate_vv: %d steps\n",this_node,n_steps);
#endif

  /* check init */
  if(rebuild_verletlist == 1) {
    sort_particles_into_cells();
    exchange_ghost();
    build_verlet_list();
  }
  if(calc_forces_first == 1) {
    force_calc();
    calc_forces_first = 0;
  }

  /* integration loop */
  for(i=0;i<n_steps;i++) {
    propagate_velocities();
    propagate_positions();
    if(rebuild_verletlist == 1) {
      sort_particles_into_cells();
      exchange_ghost();
      build_verlet_list();
    }
    else {
      exchange_ghost_pos();
    }
    force_calc();
    exchange_ghost_forces();
    propagate_velocities();
  }

}

void integrate_vv_exit()
{
#ifdef DEBUG
  if(this_node==0) fprintf(stderr,"%d: integrate_vv_exit\n",this_node);
#endif
  cells_exit();
  free(local_index);
  ghost_exit();
  verlet_exit();
  force_exit();
}

void propagate_velocities() 
{
#ifdef DEBUG
  if(this_node==0) fprintf(stderr,"%d: propagate_velocities:\n",this_node);
#endif
}

void propagate_positions() 
{
#ifdef DEBUG
  if(this_node==0) fprintf(stderr,"%d: propagate_positions:\n",this_node);
#endif
}

