#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "debug.h"
/* from these modules we modify variables: */
#include "communication.h"
#include "grid.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "integrate.h"
#include "thermostat.h"

/**********************************************
 * description of variables
 * callbacks please define where the variables
 * comes from.
 **********************************************/

/** Read-only callback for \ref fields.
    If you choose this, the variable cannot be
    changed by Tcl script code. */
int ro_callback(Tcl_Interp *interp, void *data);

/** Callback for box_l. Sets the box dimensions. */
int boxl_callback(Tcl_Interp *interp, void *_data);

/** Callback for gamma. Sets the friction coefficient gamma. */
int gamma_callback(Tcl_Interp *interp, void *_data);

/** List of all Tcl accessible global variables. If you
    want to add a new variable, ADD IT ALWAYS AT THE END.
    You should also add an \verbatim #define FIELD_*\endverbatim
    in \ref global.h.
*/
const Datafield fields[] = {
  {&nprocs,    TYPE_INT,    1, "nprocs",    ro_callback }, /* communication.c */
  {processor_grid, TYPE_INT, 3, "procgrid", pgrid_callback }, /* grid.c */
  {local_box_l, TYPE_DOUBLE, 3, "local_box_l", ro_callback }, /* global.c */
  {box_l, TYPE_DOUBLE, 3, "box_l", boxl_callback },
  {&n_total_particles, TYPE_INT, 1, "nparticles", ro_callback },
  {&n_particle_types, TYPE_INT, 1, "nptypes", ro_callback },
  {&n_interaction_types, TYPE_INT, 1, "niatypes", niatypes_callback },
  {&time_step, TYPE_DOUBLE, 1, "time_step", ro_callback }, /* integrator.c */
  {&max_cut, TYPE_DOUBLE,   1, "max_cut", ro_callback },
  {&skin, TYPE_DOUBLE,   1, "skin", ro_callback },
  {&max_range, TYPE_DOUBLE,   1, "max_range", ro_callback },
  {&friction_gamma, TYPE_DOUBLE,   1, "gamma", gamma_callback },
  { NULL, 0, 0, NULL, NULL }
};

/**********************************************
 * functions
 **********************************************/

int ro_callback(Tcl_Interp *interp, void *data)
{
  Tcl_AppendResult(interp, "variable is readonly", (char *)NULL);
  return (TCL_ERROR);
}

int boxl_callback(Tcl_Interp *interp, void *_data)
{
  double *data = _data;

  if ((data[0] < 0) || (data[1] < 0) || (data[2] < 0)) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }

  box_l[0] = data[0];
  box_l[1] = data[1];
  box_l[2] = data[2];

  changed_topology();

  return (TCL_OK);
}

int gamma_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;

  if (data < 0) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }
  friction_gamma = data;
  return (TCL_OK);
}

int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  double dbuffer[MAX_DIMENSION];
  int    ibuffer[MAX_DIMENSION];
  char   buffer[TCL_DOUBLE_SPACE + 5];
  int i, j;
  int status;

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
	    if (Tcl_GetInt(interp, argv[2 + j], &ibuffer[j]) == TCL_ERROR)
	      return (TCL_ERROR);
	    break;
	  case TYPE_DOUBLE:
	    if (Tcl_GetDouble(interp, argv[2 + j], &dbuffer[j]))
	      return (TCL_ERROR);
	    break;
	  default: ;
	  }
	}

	/* call changeproc */
	switch(fields[i].type) {
	case TYPE_INT:
	  status = fields[i].changeproc(interp, ibuffer);
	  break;
	case TYPE_DOUBLE:
	  status = fields[i].changeproc(interp, dbuffer);
	  break;
	default:
	  status = TCL_ERROR;
	}
	if (status != TCL_OK)
	  return TCL_ERROR;

	mpi_bcast_parameter(i);
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
