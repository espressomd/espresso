#include "field.h"
#include "global.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* get/set contents */
int mdvar(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv);

void field_var_init()
{
  int i;

  /* standard global variables */
  for (i = 0; fields[i].data != NULL; i++)
    if (fields[i].dimension != DIM_SCALAR)
      *((void **) fields[i].data) = NULL;
}

void field_var_resize(int new_size, int dim)
{
  void **dataptr;
  int i;

  /* standard global variables */
  for (i = 0; fields[i].data != NULL; i++) {
    if (fields[i].dimension == dim) {
      dataptr = (void **)fields[i].data;
      if (*dataptr)
	free(*dataptr);
      switch(fields[i].type) {
      case TYPE_INT:
	*(int **)dataptr = (int *)malloc(new_size*sizeof(int));
	break;
      case TYPE_DOUBLE:
	*(double **)dataptr = (double *)malloc(new_size*sizeof(double));
	break;
      default: break;
      }
    }
  }
}

void field_init(Tcl_Interp *interp)
{
  init_variables();

  Tcl_CreateCommand(interp, "mdvar", mdvar, 0, NULL);
}

int mdvar(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  char buffer[64];
  int index_start, index_end, index_max;
  int i, j;

  if (argc < 2) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     argv[0], " <variable> ?index? ?value?\"", (char *) NULL);
    return (TCL_ERROR);
  }

  for (i = 0; fields[i].data != NULL; i++) {
    if (!strncmp(argv[1], fields[i].name, strlen(argv[1]))) {      
      switch (fields[i].dimension) {
      case DIM_SCALAR:
	/* set */
	if (argc == 3) {
	  switch(fields[i].type) {
	  case TYPE_INT:
	    if (fields[i].changeproc)
	      ((IntVarChangeProc *)fields[i].changeproc)(interp, atol(argv[2]));
	    else
	      *(int *)fields[i].data = atol(argv[2]);
	    break;
	  case TYPE_DOUBLE:
	    if (fields[i].changeproc)
	      ((DoubleVarChangeProc *)fields[i].changeproc)(interp, atof(argv[2]));
	    else	  
	      *(double *)fields[i].data = atof(argv[2]);
	    break;
	  default: ;
	  }
	  return (TCL_OK);
	}

	/* get */
	switch (fields[i].type) {
	case TYPE_INT:
	  sprintf(buffer, "%d", *(int *)fields[i].data);
	  break;
	case TYPE_DOUBLE:
	  sprintf(buffer, "%10.6e", *(double *)fields[i].data);
	  break;
	default: ;
	}
	Tcl_AppendResult(interp, buffer, (char *) NULL);       
	return (TCL_OK);
      case DIM_NPARTICLES:
	index_max   = N;
	index_start = 0;
	index_end   = N;
	break;
      case DIM_NCHAINS:
	index_max   = Nchains;
	index_start = 0;
	index_end   = N;
	break;
      default:
	Tcl_AppendResult(interp, "internal error", (char *)NULL);
	return (TCL_ERROR);
      }

      /* get index range */
      if (argc >= 3) {
	switch (sscanf(argv[2], "%d-%d", &index_start, &index_end)) {
	case 1:
	  index_end = index_start + 1;
	  break;
	case 2:
	  break;
	default:
	  Tcl_AppendResult(interp, "couldn't parse index", (char *)NULL);
	  return (TCL_ERROR);
	}
      }

      /* cut indices */
      if (index_end > index_max)
	index_end = index_max;
      if (index_start >= index_end) {
	Tcl_AppendResult(interp, "index \"", argv[2], "\" out of range",
			 (char *)NULL);
	return (TCL_ERROR);
      }
      if ((index_start < 0) || (index_end < 0)) {
	Tcl_AppendResult(interp, "index < 0 ??", (char *)NULL);
	return (TCL_ERROR);
      }

      /* set */
      if (argc == 4) {
	if (index_end != index_start + 1) {
	  Tcl_AppendResult(interp, "mdvar only sets scalar values, "
			   "but range is \"", argv[2], "\"",
			   (char *)NULL);
	  return (TCL_ERROR);
	}

	switch(fields[i].type) {
	case TYPE_INT:
	  if (fields[i].changeproc)
	    ((IntVarChangeProc *)fields[i].changeproc)(interp, atol(argv[3]));
	  else
	    (*(int **)fields[i].data)[index_start] = atol(argv[3]);
	  break;
	case TYPE_DOUBLE:
	  if (fields[i].changeproc)
	    ((DoubleVarChangeProc *)fields[i].changeproc)(interp, atof(argv[3]));
	  else	  
	    (*(double **)fields[i].data)[index_start] = atof(argv[3]);
	  break;
	default: ;
	}
	return (TCL_OK);
      }

      /* get */
      for (j = index_start; j < index_end; j++) {
	switch (fields[i].type) {
	case TYPE_INT:
	  sprintf(buffer, "%d", (*(int **)fields[i].data)[j]);
	  break;
	case TYPE_DOUBLE:
	  sprintf(buffer, "%10.6e", (*(double **)fields[i].data)[j]);
	  break;
	default: ;
	}
	if (j < index_end - 1)
	  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
	else
	  Tcl_AppendResult(interp, buffer, (char *) NULL);       
      }
      return (TCL_OK);
    }
  }

  Tcl_AppendResult(interp, "variable not found \"",
		   argv[1], "\"", (char *) NULL);
  return (TCL_ERROR);
}
