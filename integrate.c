/*******************  INTEGRATE.C  *******************/
#include <stdio.h>
#include <stdlib.h>
#include "communication.h"
#include "global.h"

#include "integrate.h"
#include "cells.h"
#include "verlet.h"
#include "forces.h"

#define DEBUG

double time_step = 0.001;
double max_cut = 2.0;
double skin = 0.4;
double max_range;
double max_range2;


void tcl_integrator_init(Tcl_Interp *interp)
{
  Tcl_CreateCommand(interp, "integrate", integrate, 0, NULL);
}

int integrate(ClientData data, Tcl_Interp *interp,
	      int argc, char **argv)
{
  int  succes;
  int  n_steps;
  char buffer[256];
  
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


int integrate_vv_init()
{


#ifdef DEBUG
  if(this_node==0) fprintf(stderr,"%d: integrate_vv_init\n",this_node);
#endif
  max_range  = max_cut + skin;
  max_range2 = max_range* max_range;

  cells_init();

}

int integrate_vv(int n_steps)
{
#ifdef DEBUG
  if(this_node==0) fprintf(stderr,"%d: integrate_vv: %d steps\n",this_node,n_steps);
#endif
}

int integrate_vv_exit()
{
#ifdef DEBUG
  if(this_node==0) fprintf(stderr,"%d: integrate_vv_exit\n",this_node);
#endif
}
