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
/** tcl procedure to broadcast system parameters. */
int set_syspar(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv);

/** callback for n_total_particles */
int npart_callback(Tcl_Interp *interp, void *data);
/** callback for ro */
int ro_callback(Tcl_Interp *interp, void *data);
/** callback for procgrid */
int pgrid_callback(Tcl_Interp *interp, void *data);
/** callback for box_l */
int boxl_callback(Tcl_Interp *interp, void *data);

/** maximal length of a writable datafield */
#define MAX_DIMENSION 3
Tcl_Datafield fields[] = {
  {&nprocs,    TYPE_INT,    1, "nprocs",    ro_callback },
  {processor_grid, TYPE_INT, 3, "procgrid", pgrid_callback },
  {neighbors, TYPE_INT,    6, "neighbors", ro_callback },
  {box_l,     TYPE_DOUBLE, 3, "box_l",     boxl_callback },
  {my_left,   TYPE_DOUBLE, 3, "my_left",   ro_callback },
  {my_right,  TYPE_DOUBLE, 3, "my_right",  ro_callback },
  {&n_total_particles, TYPE_INT, 1, "nparticles", npart_callback },
  { NULL, 0, 0, NULL, NULL }
};

void tcl_datafield_init(Tcl_Interp *interp)
{
  Tcl_CreateCommand(interp, "setmd", setmd, 0, NULL);
  Tcl_CreateCommand(interp, "part", part, 0, NULL);
  Tcl_CreateCommand(interp, "set_syspar", set_syspar, 0, NULL);
}

int setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  double dbuffer[MAX_DIMENSION];
  int    ibuffer[MAX_DIMENSION];
  char   buffer[256];
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
	  return fields[i].changeproc(interp, ibuffer);
	  break;
	case TYPE_DOUBLE:
	  return fields[i].changeproc(interp, dbuffer);
	  break;
	default: return (TCL_OK);
	}
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

/* callback for npart */
int npart_callback(Tcl_Interp *interp, void *data)
{
  n_total_particles = *(int *)data;
  return (TCL_OK);
}

/* callback for ro */
int ro_callback(Tcl_Interp *interp, void *data)
{
  Tcl_AppendResult(interp, "variable is readonly", (char *)NULL);
  return (TCL_ERROR);
}

/* callback for procgrid */
int pgrid_callback(Tcl_Interp *interp, void *_data)
{
  int *data = (int *)_data;
  if ((data[0] < 0) || (data[1] < 0) || (data[2] < 0)) {
    Tcl_AppendResult(interp, "illegal value", (char *) NULL);
    return (TCL_ERROR);
  }
  processor_grid[0] = data[0];
  processor_grid[1] = data[1];
  processor_grid[2] = data[2];
  if (!setup_processor_grid()) {
    Tcl_AppendResult(interp, "processor grid does not fit nprocs", (char *) NULL);
    return (TCL_ERROR);
  }
  rebuild_verletlist = 1;
  return (TCL_OK);
}

/* callback for box_l */
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

  rebuild_verletlist = 1;

  return (TCL_OK);
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
    for (j = 0; j < part.n_bond; j++) {
      sprintf(buffer, "{%d %d}", part.bonds[j], part.bond_type[j]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      if (j < part.n_bond - 1)
	Tcl_AppendResult(interp, " ", (char *)NULL);
    }
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
      /* wird eh gefaltet! */
      //if ((pos[j] < 0) || (pos[j] >= box_l[j])) {
      //Tcl_AppendResult(interp, "particle out of bounds", (char *)NULL);
      //return (TCL_ERROR);
      //}
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

int set_syspar(ClientData data, Tcl_Interp *interp,
	       int argc, char **argv)
{
  if (argc > 1) {
    Tcl_AppendResult(interp, "No arguments expected. Any argument will be ignored ", (char *) NULL);
  }
  
  mpi_set_system_parameters();
  return (TCL_OK);
}
