/*****************************************************/
/*******************  INTEGRATE.C  *******************/
/*****************************************************/

#include "integrate.h"
#include "debug.h"

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
  if (argc == 2) {
    mpi_integrate(n_steps);
    return (TCL_OK);
  }
  else {
    Tcl_AppendResult(interp, "too many arguments:  should be \"",
		     argv[0], " <task> \"", (char *) NULL);
    return (TCL_ERROR);
  }
}


void integrate_vv_init()
{
  int i;

  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv_init\n",this_node));
  INTEG_TRACE(fprintf(stderr,"%d: nproc =%d npart=%d\n",
		      this_node,nprocs,n_total_particles));
  INTEG_TRACE(fprintf(stderr,"%d: n_total_part=%d  n_particles = %d\n",
		      this_node,n_total_particles,n_particles));
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
  int i;
  INTEG_TRACE(fprintf(stderr,"%d: integrate_vv: %d steps\n",this_node,
		      n_steps));

  MPI_Barrier(MPI_COMM_WORLD);

  /* check init */
  if(rebuild_verletlist == 1) {
    exchange_part();
    MPI_Barrier(MPI_COMM_WORLD);
    sort_particles_into_cells();
    MPI_Barrier(MPI_COMM_WORLD);
    exchange_ghost();
    MPI_Barrier(MPI_COMM_WORLD);
    build_verlet_list();
  }
  if(calc_forces_first == 1) {
    force_calc();
    collect_ghost_forces();
    calc_forces_first = 0;
  }

  /* integration loop */
  for(i=0;i<n_steps;i++) {
    INTEG_TRACE(fprintf(stderr,"%d: STEP %d\n",this_node,i));
    propagate_velocities();
    propagate_positions();
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

void propagate_velocities() 
{
  INTEG_TRACE(fprintf(stderr,"%d: propagate_velocities:\n",this_node));
}

void propagate_positions() 
{
  INTEG_TRACE(fprintf(stderr,"%d: propagate_positions:\n",this_node));
}

