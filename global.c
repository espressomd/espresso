#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "debug.h"
#include "communication.h"
#include "grid.h"

/**********************************************
 * description of variables
 * callbacks please define where the variables
 * comes from.
 **********************************************/

/** read only callback. If you choose this, the
    variable cannot be changed from Tcl */
int ro_callback(Tcl_Interp *interp, void *data);

/* callback for box_l */
int boxl_callback(Tcl_Interp *interp, void *_data);

/* do not change order !!!
 * and if you add something, also add an #define FIELD_* in
 * global.h
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
  { NULL, 0, 0, NULL, NULL }
};

/**********************************************
 * variables
 **********************************************/

/* simulation box and domain decompostion */ 
double box_l[3]       = {1, 1, 1};
double local_box_l[3] = {1, 1, 1};
double my_left[3]     = {0, 0, 0};
double my_right[3]    = {1, 1, 1};

/* particles */
int n_total_particles = 0;
int *particle_node = NULL;

int     n_particles = 0;
int   max_particles = 0;
int        n_ghosts = 0;
Particle *particles = NULL;

int *local_index;

/* nonbonded (short range) interactions */
int n_particle_types = 0;
int n_interaction_types = 0;
IA_parameters *ia_params = NULL;

/**********************************************
 * procedures
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
