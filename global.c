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

/**********************************************
 * description of variables
 * callbacks please define where the variables
 * comes from.
 **********************************************/

/** Read-only callback for \ref #fields.
    If you choose this, the variable cannot be
    changed by Tcl script code. */
int ro_callback(Tcl_Interp *interp, void *data);

/** List of all Tcl accessible global variables. If you
    want to add a new variable, ADD IT ALWAYS AT THE END.
    You should also add an \verbatim #define FIELD_*\endverbatim
    in \ref global.h and a descriptive text in \ref variables_page.
*/
const Datafield fields[] = {
  {box_l,            TYPE_DOUBLE, 3, "box_l",         boxl_callback, 1 },            /* grid.c */
  {cell_grid,           TYPE_INT, 3, "cell_grid",     ro_callback, 6 },              /* cells.c */
  {cell_size,        TYPE_DOUBLE, 3, "cell_size",     ro_callback, 6 },              /* cells.c */
  {&friction_gamma,  TYPE_DOUBLE, 1, "gamma",         gamma_callback, 1 },           /* thermostat.c */
  {&lj_force_cap,    TYPE_DOUBLE, 1, "lj_force_cap",  lj_force_cap_callback, 2 },    /* interaction.c */
  {local_box_l,      TYPE_DOUBLE, 3, "local_box_l",   ro_callback, 2 },              /* global.c */
  {&max_cut,         TYPE_DOUBLE, 1, "max_cut",       ro_callback, 5 },              /* interaction_data.c */
  {&max_num_cells,      TYPE_INT, 1, "max_num_cells", max_num_cells_callback, 5 },   /* cells.c */
  {&max_seen_particle,  TYPE_INT, 1, "max_part",      ro_callback, 5 },              /* particle_data.c */
  {&max_range,       TYPE_DOUBLE, 1, "max_range",     ro_callback, 5 },              /* integrate.c */
  {&max_skin,        TYPE_DOUBLE, 1, "max_skin",      ro_callback, 5 },              /* integrate.c */
  {&n_nodes,            TYPE_INT, 1, "n_nodes",       ro_callback, 3 },              /* communication.c */
  {&n_total_particles,  TYPE_INT, 1, "n_part",        ro_callback, 6 },              /* particle.c */
  {&n_particle_types,   TYPE_INT, 1, "n_part_types",  ro_callback, 8 },              /* interaction_data.c */
  {node_grid,           TYPE_INT, 3, "node_grid",     node_grid_callback, 2 },       /* grid.c */
#ifdef PARTIAL_PERIODIC
  {periodic,            TYPE_INT, 3, "periodicity",   per_callback, 1 },             /* grid,c */
#else
  {periodic,            TYPE_INT, 3, "periodicity",   ro_callback, 1 },              /* grid,c */
#endif
  {&skin,            TYPE_DOUBLE, 1, "skin",          skin_callback, 2 },            /* integrate.c */
  {&temperature,     TYPE_DOUBLE, 1, "temperature",   temp_callback, 2 },            /* thermostat.c */
  {&sim_time,        TYPE_DOUBLE, 1, "time",          start_time_callback, 4 },      /* integrate.c */
  {&time_step,       TYPE_DOUBLE, 1, "time_step",     time_step_callback, 5 },       /* integrate.c */
  {&transfer_rate,      TYPE_INT, 1, "transfer_rate", ro_callback, 2 }     ,         /* imd.c */
  {&rebuild_verletlist, TYPE_INT, 1, "verlet_flag",   rebuild_vlist_callback, 8 },   /* verlet.c */
  {&verlet_reuse,    TYPE_DOUBLE, 1, "verlet_reuse",  ro_callback, 8 },              /* integrate.c */
  { NULL, 0, 0, NULL, NULL, 0 }
};

/** \page variables_page Global variables

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
	     \ref temperature - Temperature of the simulation.
	     Enters the thermostat and the coulomb prefactor = bjerrum * temperature.
	<li> \verbatim time double \endverbatim
	     \ref time - The simulation time.
	<li> \verbatim time_step double \endverbatim
	     \ref time_step - Time step for MD integration
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
  Tcl_AppendResult(interp, "Warning: variable is readonly", (char *)NULL);
  return (TCL_OK);
}

int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  char databuf[MAX_DIMENSION*(sizeof(int) + sizeof(double))];
  char buffer[TCL_DOUBLE_SPACE + 5];
  int i, j;

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <variable> ?value? ?value? ...\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  for (i = 0; fields[i].data != NULL; i++) {
    if (!strncmp(argv[1], fields[i].name, strlen(argv[1]))) {
      if(strlen(argv[1]) < fields[i].min_length) {
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

      /* get */
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

      return (TCL_OK);
    }
  }
  Tcl_AppendResult(interp, "unknown variable \"",
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
