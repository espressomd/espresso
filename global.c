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

/** Read-only callback for \ref fields.
    If you choose this, the variable cannot be
    changed by Tcl script code. */
int ro_callback(Tcl_Interp *interp, void *data);

/** List of all Tcl accessible global variables. If you
    want to add a new variable, ADD IT ALWAYS AT THE END.
    You should also add an \verbatim #define FIELD_*\endverbatim
    in \ref global.h and a descriptive text in \ref variables_page.
*/
const Datafield fields[] = {
  {&n_nodes,    TYPE_INT,    1, "n_nodes",    ro_callback },           /* communication.c */
  {node_grid, TYPE_INT, 3, "node_grid", node_grid_callback },          /* grid.c */
  {local_box_l, TYPE_DOUBLE, 3, "local_box_l", ro_callback },          /* global.c */
  {box_l, TYPE_DOUBLE, 3, "box_l", boxl_callback },                    /* grid.c */
  {&max_seen_particle, TYPE_INT, 1, "maxpart", ro_callback },          /* particle_data.c */
  {&n_particle_types, TYPE_INT, 1, "nptypes", ro_callback },           /* interaction_data.c */
  {&time_step, TYPE_DOUBLE, 1, "time_step", time_step_callback },      /* integrate.c */
  {&max_cut, TYPE_DOUBLE,   1, "max_cut", ro_callback },               /* interaction_data.c */
  {&skin, TYPE_DOUBLE,   1, "skin", skin_callback },                   /* integrate.c */
  {&max_range, TYPE_DOUBLE,   1, "max_range", ro_callback },           /* integrate.c */
  {&friction_gamma, TYPE_DOUBLE,   1, "gamma", gamma_callback },       /* thermostat.c */
  {&rebuild_verletlist, TYPE_INT, 1, "verletflag", rebuild_vlist_callback }, /* verlet.c */
  {&(p3m.bjerrum), TYPE_DOUBLE,   1, "bjerrum", bjerrum_callback },    /* p3m.c */
  {&(p3m.alpha), TYPE_DOUBLE,   1, "p3m_alpha", p3malpha_callback },   /* p3m.c */
  {&(p3m.r_cut), TYPE_DOUBLE,   1, "p3m_r_cut", p3mrcut_callback },    /* p3m.c */
  {p3m.mesh, TYPE_INT,   3, "p3m_mesh", p3mmesh_callback },            /* p3m.c */
  {&(p3m.cao), TYPE_INT,   2, "p3m_cao", p3mcao_callback },            /* p3m.c */
  {&(p3m.epsilon), TYPE_DOUBLE,   1, "p3m_epsilon", p3mepsilon_callback },    /* p3m.c */
  {p3m.mesh_off, TYPE_DOUBLE,   3, "p3m_mesh_offset", p3mmeshoff_callback },  /* p3m.c */
  {&transfer_rate, TYPE_INT,   1, "transfer_rate", ro_callback },             /* imd.c */
  {&max_num_cells, TYPE_INT,   1, "max_num_cells", max_num_cells_callback },  /* cells.c */
#ifdef PARTIAL_PERIODIC
  {periodic, TYPE_INT,   3, "periodicity", per_callback },             /* grid,c */
#else
  {periodic, TYPE_INT,   3, "periodicity", ro_callback },              /* grid,c */
#endif
  {&temperature, TYPE_DOUBLE, 1, "temp", temp_callback },              /* thermostat.c */
  {&lj_force_cap, TYPE_DOUBLE, 1, "lj_force_cap", lj_force_cap_callback },  /* interaction.c */
  {&start_time, TYPE_DOUBLE, 1, "start_time", start_time_callback }, /* integrate.c */
  {&sim_time, TYPE_DOUBLE, 1, "time", ro_callback },                 /* integrate.c */
  { NULL, 0, 0, NULL, NULL }
};

/** \page variables_page Global variables

     The following list explains the usage of the variables that are
     accessible via \ref tcl_setmd.  The list gives the setmd name of
     (hopefully) all available variables, the data type and a link to
     the documentation of the corresponding C variable.
     Variables that are marked read only can only be written by C code.

	<ul>
	<li> \verbatim n_nodes int (ro) \endverbatim
	\ref n_nodes - Number of nodes.
	<li> \verbatim node_grid int[3] \endverbatim
	\ref node_grid - 3D node grid for real space domain
	decomposition (optional,
	if unset an optimal set is chosen automatically).	
	<li> \verbatim local_box_l int[3] (ro) \endverbatim
	\ref local_box_l - Local simulation box length of the nodes.
	<li> \verbatim box_l double[3] \endverbatim
	\ref box_l - Simulation box length.
	<li> \verbatim maxpart int (ro) \endverbatim
	\ref max_seen_particle - Maximal identity of a particle.
	THIS IS IN GENERAL _NOT_ RELATED
	TO THE NUMBER OF PARTICLES.
	<li> \verbatim nptypes int (ro) \endverbatim
	\ref n_particle_types - Number of particle
	types that were used so far in \ref tcl_inter.
	<li> \verbatim niatypes int \endverbatim
	\ref n_interaction_types - Number of interaction types
	(not fully implemented. DO NOT USE).
	<li> \verbatim time_step double (_currently_ ro)\endverbatim
	\ref time_step - Time step for MD integration
	(_currently_ read only, set to 0.001).
	<li> \verbatim max_cut double (ro) \endverbatim
	\ref max_cut - Maximal cutoff of real space interactions.
	<li> \verbatim skin double \endverbatim
	\ref skin - Skin for the Verlet list.
	<li> \verbatim max_range double (ro)\endverbatim
	\ref max_range - Maximal range of real space interactions: max_cut + skin.
	<li> \verbatim gamma double \endverbatim
	\ref friction_gamma - Friction constant.
	<li> \verbatim verletflag bool (ro)\endverbatim
	\ref rebuild_verletlist - Indicates wether the Verlet list will be rebuild.
	The program decides this automatically based on your actions on the data.
	<li> \verbatim bjerrum double \endverbatim
	\ref p3m_struct::bjerrum - Bjerrum length. If 0, electrostatic interaction is turned off.
	<li> \verbatim p3m_alpha double \endverbatim
	\ref p3m_struct::alpha - Ewald splitting parameter.
	<li> \verbatim p3m_r_cut double \endverbatim
	\ref p3m_struct::r_cut - Real space cutoff for P3M.
	<li> \verbatim p3m_mesh int[3]\endverbatim
	\ref p3m_struct::mesh - Mesh size for P3M k space.
	<li> \verbatim p3m_cao int \endverbatim
	\ref p3m_struct::cao - Charge assignment order for P3M particle-mesh interaction (1 < 
	p3m_cao < 7).
	<li> \verbatim p3m_epsilon double\endverbatim
	\ref p3m_struct::epsilon - Dielectric constant at infinity
	(boundary condition for Ewald sumation).
	<li> \verbatim p3m_mesh_offset double[3] \endverbatim
	\ref p3m_struct::mesh_off - Offset of the first mesh point from the origin in mesh
	coordinates (all values between 0.0 and 1.0).
	<li> \verbatim transfer_rate int (ro)\endverbatim
	\ref transfer_rate - Tranfer rate for VMD connection. You can use this
	to transfer any integer value to the simulation from VMD.
	<li> \verbatim max_num_cells int> \endverbatim
	\ref max_num_cells - Maximal number of cells for the link cell
	algorithm. Reasonable values are between 125 and 1000, or for
	some problems (n_total_particles/n_nodes).
	<li> \verbatim periodicity bool[3]\endverbatim
	\ref periodic - Specifies periodicity for the three directions.
	This variable is read-only and returns (1,1,1) without the compiler flag
	\ref PARTIAL_PERIODIC from \ref config.h "config.h" .
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

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <variable> ?value? ?value? ...\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  for (i = 0; fields[i].data != NULL; i++) {
    if (!strncmp(argv[1], fields[i].name, strlen(argv[1]))) {
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
  Tcl_AppendResult(interp, "unkown variable \"",
		   argv[1], "\"", (char *) NULL);
  return (TCL_ERROR);
}
