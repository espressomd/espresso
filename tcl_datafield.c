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
/** tcl procedure for interaction access */
int inter(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);
/** tcl procedure for particle access */
int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

void tcl_datafield_init(Tcl_Interp *interp)
{
  Tcl_CreateCommand(interp, "setmd", setmd, 0, NULL);
  Tcl_CreateCommand(interp, "inter", inter, 0, NULL);
  Tcl_CreateCommand(interp, "part", part, 0, NULL);
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
  if (part_num < 0) {
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
    /* FIXME: new structure
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
  
  /* set particle data */
  argc -= 2;
  argv += 2;
  while (argc > 0) {
    if (!strncmp(argv[0], "pos", strlen(argv[0]))) {
      double pos[3];
      if (argc < 4) {
	Tcl_AppendResult(interp, "pos requires 3 arguments", (char *) NULL);
	return (TCL_ERROR);
      }
      /* set position */
      for (j = 0; j < 3; j++) {
	if (Tcl_GetDouble(interp, argv[1 + j], &pos[j]) == TCL_ERROR)
	  return (TCL_ERROR);
      }
    
      if (node == -1) {
	/* spatial decomposite position */
	node = find_node(pos);
	mpi_attach_particle(part_num, node);
      }

      mpi_send_pos(node, part_num, pos);

      argc -= 4;
      argv += 4;
    }
    else {
      if (node == -1) {
	Tcl_AppendResult(interp, "set particle position first", (char *)NULL);
	return (TCL_ERROR);
      }
      
      if (!strncmp(argv[0], "q", strlen(argv[0]))) {
	double q;
	if (argc < 2) {
	  Tcl_AppendResult(interp, "q requires 1 argument", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set charge */
	if (Tcl_GetDouble(interp, argv[1], &q) == TCL_ERROR)
	  return (TCL_ERROR);
	
	mpi_send_q(node, part_num, q);

	argc -= 2;
	argv += 2;
      }
      else if (!strncmp(argv[0], "type", strlen(argv[0]))) {
	int type;
	if (argc < 2) {
	  Tcl_AppendResult(interp, "type requires 1 argument", (char *) NULL);
	  return (TCL_ERROR);
	}
	/* set type */
	if (Tcl_GetInt(interp, argv[1], &type) == TCL_ERROR)
	  return (TCL_ERROR);

	if (type < 0) {
	  Tcl_AppendResult(interp, "invalid particle type", (char *) NULL);
	  return (TCL_ERROR);	  
	} 

	// make sure type exists
	realloc_ia_params(type);

	mpi_send_type(node, part_num, type);

	argc -= 2;
	argv += 2;
      }
      else {
	Tcl_AppendResult(interp, "unknown particle parameter \"", argv[0],"\"", (char *)NULL);
	return (TCL_ERROR);
      }
    }
  }

  return (TCL_OK);
}

int inter(ClientData _data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  int i, j;
  IA_parameters *data, *data_sym;

  if (argc < 3) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <type 1> <type 2> ?interaction? ?values?\"",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if ((Tcl_GetInt(interp, argv[1], &i) == TCL_ERROR) ||
      (Tcl_GetInt(interp, argv[2], &j) == TCL_ERROR))
    return (TCL_ERROR);

  data     = safe_get_ia_param(i, j);
  data_sym = safe_get_ia_param(j, i);

  if (!data || !data_sym) {
    Tcl_AppendResult(interp, "particle types must be nonnegative",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  if (argc == 3) {
    /* print interaction information */
    char buffer[256];
    sprintf(buffer, "{lennard-jones %10.5e %10.5e %10.5e %10.5e}",
	    data->LJ_epsilon, data->LJ_cutoff,
	    data->LJ_shift, data->LJ_offset);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
  }

  /* set interaction parameters */
  argc -= 3;
  argv += 3;

  while (argc > 0) {
    if (!strncmp(argv[0], "lennard-jones", strlen(argv[0]))) {
      if (argc < 5) {
	Tcl_AppendResult(interp, "lennard-jones needs 4 parameters: "
			 "<lj_epsilon> <lj_cut> <lj_shift> <lj_offset>",
			 (char *) NULL);
	return (TCL_ERROR);
      }
      
      if ((Tcl_GetDouble(interp, argv[1], &data->LJ_epsilon) == TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[2], &data->LJ_cutoff)  == TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[3], &data->LJ_shift)   == TCL_ERROR) ||
	  (Tcl_GetDouble(interp, argv[4], &data->LJ_offset)  == TCL_ERROR))
	return (TCL_ERROR);

      /* LJ should be symmetrically */
      data_sym->LJ_epsilon = data->LJ_epsilon;
      data_sym->LJ_cutoff  = data->LJ_cutoff;
      data_sym->LJ_shift   = data->LJ_shift;
      data_sym->LJ_offset  = data->LJ_offset;
      argc -= 5;
      argv += 5;

      mpi_bcast_ia_params(i, j);
    }
    else {
      Tcl_AppendResult(interp, "unknown interaction type \"", argv[3],
		       "\"", (char *)NULL);
      return (TCL_ERROR);
    }
  }

  return (TCL_OK);
}
