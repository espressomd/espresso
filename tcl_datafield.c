#include "tcl_datafield.h"
#include "global.h"
#include "communication.h"
#include "grid.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* cwz-build-comman: ssh chakotay "builtin cd /nhomes/janeway/axel/progs/tcl_md; make" 
   cwz-build-command: make
*/

/** tcl procedure for datafield access */
int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);
/** tcl procedure for particle access */
int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);
/* tcl procedure to broadcast system parameters.
int set_syspar(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv);
*/

void tcl_datafield_init(Tcl_Interp *interp)
{
  Tcl_CreateCommand(interp, "setmd", setmd, 0, NULL);
  Tcl_CreateCommand(interp, "part", part, 0, NULL);
  //  Tcl_CreateCommand(interp, "set_syspar", set_syspar, 0, NULL);
}

int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  double dbuffer[MAX_DIMENSION];
  int    ibuffer[MAX_DIMENSION];
  char   buffer[256];
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
	  sprintf(buffer, "%10.6e", ((double *)fields[i].data)[j]);
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

int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv)
{
  int part_num = -1;
  int node, j;
  char buffer[256];

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <part num> ?what? ?value?\"", (char *) NULL);
    return (TCL_ERROR);
  }
  if (!processor_grid_is_set())
    setup_processor_grid();

  part_num = atol(argv[1]);
  if ((part_num < 0) || (part_num >= n_total_particles)) {
    Tcl_AppendResult(interp, "illegal particle", (char *) NULL);
    return (TCL_ERROR);
  }

  node = mpi_who_has(part_num);

  /* print out particle information */
  if (argc == 2) {
    Particle part;
    /* retrieve particle data */
    if (node == -1) {
      Tcl_AppendResult(interp, "particle not found", (char *) NULL);
      return (TCL_ERROR);
    }
    mpi_recv_part(node, part_num, &part);
    sprintf(buffer, "%d %d {%7.3e %7.3e %7.3e} %8.2e \n",
	    node, part.type,
	    part.p[0], part.p[1], part.p[2],
	    part.q);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    sprintf(buffer, "       {%7.3e %7.3e %7.3e} {%7.3e %7.3e %7.3e}",
	    part.v[0], part.v[1], part.v[2],
	    part.f[0], part.f[1], part.f[2]);
    Tcl_AppendResult(interp, buffer, " {", (char *)NULL);
    /*
    for (j = 0; j < part.n_bonds; j++) {
      sprintf(buffer, "{%d %d}", part.bonds[j], part.bond_type[j]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      if (j < part.n_bond - 1)
	Tcl_AppendResult(interp, " ", (char *)NULL);
    }
    */
    Tcl_AppendResult(interp, "}", (char *)NULL);
    return (TCL_OK);
  }
  
  if (!strncmp(argv[2], "pos", strlen(argv[2]))) {
    double pos[3];
    if (argc != 6) {
      Tcl_AppendResult(interp, "pos requires 3 arguments");
      return (TCL_ERROR);
    }
    /* set position */
    for (j = 0; j < 3; j++) {
      if (Tcl_GetDouble(interp, argv[3 + j], &pos[j]) == TCL_ERROR)
	return (TCL_ERROR);
      pos[j] -= floor(pos[j]/box_l[j])*box_l[j];
    }
    
    if (node == -1) {
      /* spatial decomposite position */
      node = find_node(pos);
      mpi_attach_particle(part_num, node);
    }

    mpi_send_pos(node, part_num, pos);
    return (TCL_OK);
  }

  if (node == -1) {
    Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
    return (TCL_ERROR);
  }

  Tcl_AppendResult(interp, "unknown job \"", argv[2],"\"", (char *)NULL);
  return (TCL_ERROR);
}

/*
int set_syspar(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv)
{
  if (argc > 1) {
    Tcl_AppendResult(interp, "No arguments expected. Any argument will be ignored ", (char *) NULL);
  }
  
  mpi_set_system_parameters();
  return (TCL_OK);
}
*/
