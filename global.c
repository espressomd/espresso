/** \file global.c
    Implementation of \ref global.h "global.h".
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "global.h"
#include "debug.h"
/* from these modules we modify variables: */
#include "communication.h"
#include "cells.h"
#include "grid.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "integrate.h"
#include "thermostat.h"
#include "forces.h"
#include "verlet.h"
#include "p3m.h"
#include "imd.h"
#include "tuning.h"

/**********************************************
 * description of variables
 * callbacks please define where the variables
 * comes from.
 **********************************************/

/** Read-only callback for \ref #fields.
    If you choose this, the variable cannot be
    changed by Tcl script code. */
int ro_callback(Tcl_Interp *interp, void *data);

/// List of all Tcl accessible global variables
const Datafield fields[] = {
  {box_l,            TYPE_DOUBLE, 3, "box_l",         boxl_callback, 1 },            /* 0 from grid.c */
  {cell_grid,           TYPE_INT, 3, "cell_grid",     ro_callback, 6 },              /* 1 from cells.c */
  {cell_size,        TYPE_DOUBLE, 3, "cell_size",     ro_callback, 6 },              /* 2 from cells.c */
  {&friction_gamma,  TYPE_DOUBLE, 1, "gamma",         gamma_callback, 1 },           /* 3 from thermostat.c */
  {local_box_l,      TYPE_DOUBLE, 3, "local_box_l",   ro_callback, 2 },              /* 4 from global.c */
  {&max_cut,         TYPE_DOUBLE, 1, "max_cut",       ro_callback, 5 },              /* 5 from interaction_data.c */
  {&max_num_cells,      TYPE_INT, 1, "max_num_cells", max_num_cells_callback, 5 },   /* 6 from cells.c */
  {&max_seen_particle,  TYPE_INT, 1, "max_part",      ro_callback, 5 },              /* 7 from particle_data.c */
  {&max_range,       TYPE_DOUBLE, 1, "max_range",     ro_callback, 5 },              /* 8 from integrate.c */
  {&max_skin,        TYPE_DOUBLE, 1, "max_skin",      ro_callback, 5 },              /* 9 from integrate.c */
  {&n_nodes,            TYPE_INT, 1, "n_nodes",       ro_callback, 3 },              /* 10 from communication.c */
  {&n_total_particles,  TYPE_INT, 1, "n_part",        ro_callback, 6 },              /* 11 from particle.c */
  {&n_particle_types,   TYPE_INT, 1, "n_part_types",  ro_callback, 8 },              /* 12 from interaction_data.c */
  {node_grid,           TYPE_INT, 3, "node_grid",     node_grid_callback, 2 },       /* 13 from grid.c */
#ifdef PARTIAL_PERIODIC
  {periodic,            TYPE_INT, 3, "periodicity",   per_callback, 1 },             /* 14 from grid.c */
#else
  {periodic,            TYPE_INT, 3, "periodicity",   ro_callback, 1 },              /* 14 from grid,c */
#endif
  {&skin,            TYPE_DOUBLE, 1, "skin",          skin_callback, 2 },            /* 15 from integrate.c */
  {&temperature,     TYPE_DOUBLE, 1, "temperature",   temp_callback, 2 },            /* 16 from thermostat.c */
  {&sim_time,        TYPE_DOUBLE, 1, "time",          time_callback, 4 },            /* 17 from integrate.c */
  {&time_step,       TYPE_DOUBLE, 1, "time_step",     time_step_callback, 5 },       /* 18 from integrate.c */
  {&timing_samples,     TYPE_INT, 1, "timings",       timings_callback, 4 },         /* 19 from tuning.c */
  {&transfer_rate,      TYPE_INT, 1, "transfer_rate", ro_callback, 2 }     ,         /* 20 from imd.c */
  {&rebuild_verletlist, TYPE_INT, 1, "verlet_flag",   ro_callback, 8 },              /* 21 from verlet.c */
  {&verlet_reuse,    TYPE_DOUBLE, 1, "verlet_reuse",  ro_callback, 8 },              /* 22 from integrate.c */
  { NULL, 0, 0, NULL, NULL, 0 }
};

/** \page variables_page Global variables
    \section desc_var Description of the global variables
     The following list explains the usage of the variables that are
     accessible via \ref tcl_setmd.  The list gives the setmd name of
     (hopefully) all available variables, the data type and a link to
     the documentation of the corresponding C variable.
     Variables that are marked read only can only be written by C code.

	<ul>
	<li> \verbatim box_l double[3] \endverbatim
             \ref box_l - Simulation box length.
	<li> \verbatim cell_grid int[3] (ro) \endverbatim
             \ref #cell_grid - dimension of the inner cell grid.
	<li> \verbatim cell_size double[3] (ro) \endverbatim
	     \ref #cell_size - box length of a cell.
	<li> \verbatim gamma double \endverbatim
	     \ref friction_gamma - Friction constant.

	<li> \verbatim local_box_l int[3] (ro) \endverbatim
	     \ref local_box_l - Local simulation box length of the nodes.
	<li> \verbatim max_cut double (ro) \endverbatim
	     \ref max_cut - Maximal cutoff of real space interactions.
	<li> \verbatim max_num_cells int> \endverbatim
             \ref max_num_cells - Maximal number of cells for the link cell
	     algorithm. Reasonable values are between 125 and 1000, or for
	     some problems (\ref n_total_particles / \ref n_nodes).
	<li> \verbatim max_part int (ro) \endverbatim
  	     \ref max_seen_particle - Maximal identity of a particle.
	     THIS IS IN GENERAL _NOT_ RELATED
	     TO THE NUMBER OF PARTICLES.
	<li> \verbatim max_range double (ro)\endverbatim
	     \ref max_range - Maximal range of real space interactions: max_cut + skin.
	<li> \verbatim max_skin double (ro)\endverbatim
	     \ref max_skin - Maximal skin to be used for the link cell/verlet algorithm.
	     This is Min(\ref #cell_size) - \ref max_range.
	<li> \verbatim n_nodes int (ro) \endverbatim
  	     \ref n_nodes - Number of nodes.
	<li> \verbatim n_part int (ro) \endverbatim
  	     \ref n_total_particles - Total number of particles.
	<li> \verbatim n_part_types int (ro) \endverbatim
	     \ref n_particle_types - Number of particle
	     types that were used so far in \ref tcl_inter.
	<li> \verbatim node_grid int[3] \endverbatim
 	     \ref node_grid - 3D node grid for real space domain
	     decomposition (optional,
	     if unset an optimal set is chosen automatically).	
	<li> \verbatim periodicity bool[3]\endverbatim
             \ref #periodic - Specifies periodicity for the three directions.
             If the compiler flag PARTIAL_PERIODIC from \ref config.h "config.h" is set,
	     this variable can be set to (1,1,1) or (0,0,0) at the moment.
	     If not it is readonly and gives the default setting (1,1,1).
	<li> \verbatim skin double \endverbatim
	     \ref #skin - Skin for the Verlet list.
	<li> \verbatim temperature double \endverbatim
	     \ref #temperature - Temperature of the simulation.
	     Enters the thermostat and the coulomb prefactor = bjerrum * temperature.
	<li> \verbatim time double \endverbatim
	     \ref #time - The simulation time.
	<li> \verbatim time_step double \endverbatim
	     \ref time_step - Time step for MD integration
	<li> \verbatim timings int \endverbatim
	     \ref timing_samples - number of timing samples to take into account if set.
	<li> \verbatim transfer_rate int (ro)\endverbatim
 	     \ref transfer_rate - Transfer rate for VMD connection. You can use this
	     to transfer any integer value to the simulation from VMD.
	<li> \verbatim verlet_flag bool \endverbatim
	     \ref rebuild_verletlist - Indicates whether the Verlet list will be rebuild.
	     The program decides this normally automatically based on your actions on the data.
	     see \ref initialize.h "initialize.h" for more information
	<li> \verbatim verlet_reuse bool \endverbatim
	     \ref verlet_reuse - Average number of integration steps the verlet list has been re-used.
	</ul>    

 */

/**********************************************
 * functions
 **********************************************/


int ro_callback(Tcl_Interp *interp, void *data)
{
  Tcl_AppendResult(interp, "variable is readonly", (char *)NULL);
  return (TCL_ERROR);
}

int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  char databuf[MAX_DIMENSION*(sizeof(int) + sizeof(double))];
  char buffer[TCL_DOUBLE_SPACE + 5];
  int i, j;
  int all = (argc == 1);

  /* loop over all global variables. Has two purposes:
     either we write al variables or search for the one
     to write */
  for (i = 0; fields[i].data != NULL; i++) {
    if (all || !strncmp(argv[1], fields[i].name, strlen(argv[1]))) {
      if (!all) {
	if ((int)strlen(argv[1]) < fields[i].min_length) {
	  Tcl_AppendResult(interp, "Argument \"",argv[1],"\" not long ", (char *) NULL);
	  Tcl_AppendResult(interp, "enough to identify a setmd variable!", (char *) NULL);
	  return (TCL_ERROR);
	}
	if (argc >= 3) {
	  /* set */
	  /* parse in data */
	  if (argc != 2 + fields[i].dimension) {
	    sprintf(buffer, "%d", fields[i].dimension);	  
	    Tcl_AppendResult(interp, "\"", argv[1],
			     "\" has dimension ",
			     buffer, (char *) NULL);
	    sprintf(buffer, " not %d", argc - 2);	  
	    Tcl_AppendResult(interp, buffer, (char *) NULL);
	    return (TCL_ERROR);
	  }

	  /* get new value */
	  for (j = 0; j < fields[i].dimension; j++) {
	    switch (fields[i].type) {
	    case TYPE_INT:
	      if (Tcl_GetInt(interp, argv[2 + j], (int *)databuf + j) == TCL_ERROR)
		return (TCL_ERROR);
	      break;
	    case TYPE_DOUBLE:
	      if (Tcl_GetDouble(interp, argv[2 + j], (double *)databuf + j))
		return (TCL_ERROR);
	      break;
	    default: ;
	    }
	  }

	  if (fields[i].changeproc(interp, databuf) != TCL_OK)
	    return TCL_ERROR;
	}
      }

      /* get */
      if (all) {
	if (i != 0)
	  Tcl_AppendResult(interp, " ", (char *)NULL);
	Tcl_AppendResult(interp, "{", fields[i].name, " ", (char *)NULL);
      }
      for (j = 0; j < fields[i].dimension; j++) {
	switch (fields[i].type) {
	case TYPE_INT:
	  sprintf(buffer, "%d", ((int *)fields[i].data)[j]);
	  break;
	case TYPE_DOUBLE:
	  Tcl_PrintDouble(interp, ((double *)fields[i].data)[j], buffer);
	  break;
	default: ;
	}
	Tcl_AppendResult(interp, buffer, (char *) NULL);
	if (j < fields[i].dimension - 1)
	  Tcl_AppendResult(interp, " ", (char *) NULL);
      }
      
      if (all)
	Tcl_AppendResult(interp, "}", (char *)NULL);
      /* wrote out one value, so skip rest */
      if (!all)
	return (TCL_OK);
    }
  }

  if (all)
    return TCL_OK;

  Tcl_AppendResult(interp, "unknown md variable \"",
		   argv[1], "\"", (char *) NULL);
  return (TCL_ERROR);
}

int code_info(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  if (argc < 2) {
    version_callback(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    compilation_callback(interp);
    Tcl_AppendResult(interp, "\n", (char *) NULL);
    debug_callback(interp);
  }
  else {
    if(!strncmp(argv[1], "version" , strlen(argv[1]))) {
      version_callback(interp);
    }
    else if(!strncmp(argv[1], "compilation" , strlen(argv[1]))) {
      compilation_callback(interp);
    }
    else if(!strncmp(argv[1], "debug" , strlen(argv[1]))) {
      debug_callback(interp);
    }
    else {
      Tcl_AppendResult(interp, "info ",argv[1]," not known!", (char *) NULL);
    }
  }
  return (TCL_OK);
}
