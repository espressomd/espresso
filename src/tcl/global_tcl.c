/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file global_tcl.c
    Implementation of \ref global_tcl.h "global_tcl.h".
*/
#include "utils.h"
#include "parser.h"
#include "global.h"
#include "global_tcl.h"
#include "pressure_tcl.h"
#include "integrate_tcl.h"
#include "grid_tcl.h"
#include "domain_decomposition_tcl.h"
#include "thermostat_tcl.h"

/**********************************************
 * variables
 **********************************************/

static SetCallback **callbacks = NULL;
static int n_callbacks = 0;

/**********************************************
 * functions
 **********************************************/

/** Read-only callback for \ref #fields, the default.
    Variables without special callback cannot be written from Tcl*/
static int tclcallback_ro(Tcl_Interp *interp, void *data)
{
  Tcl_AppendResult(interp, "variable is readonly", (char *)NULL);
  return (TCL_ERROR);
}

static SetCallback *find_callback(int field)
{
  if (field >= n_callbacks)
    return tclcallback_ro;
  return callbacks[field];
}

void register_global_callback(int field, SetCallback *callback)
{
  /* resize, if necessary, and initialize fields inbetween as read-only */
  if (n_callbacks <= field) {
    callbacks = (SetCallback **)realloc(callbacks, (field + 1)*sizeof(SetCallback *));
    for (int f = n_callbacks; f < field; ++f)
      callbacks[f] = tclcallback_ro;
    n_callbacks = field + 1;
  }
  callbacks[field] = callback;
}

int tclcommand_setmd(ClientData data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  union {
    int    intbuf[MAX_DIMENSION];
    double doublebuf[MAX_DIMENSION];
  } databuf;
  char buffer[TCL_DOUBLE_SPACE + 5];
  int i, j;
  int all = (argc == 1), writing = (argc >= 3);

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
	if (writing) {
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
	      if (Tcl_GetInt(interp, argv[2 + j], databuf.intbuf + j) == TCL_ERROR)
		return (TCL_ERROR);
	      break;
	    case TYPE_BOOL: {
	      int dta;
	      if (Tcl_GetInt(interp, argv[2 + j], &dta))
		return (TCL_ERROR);
	      if (dta) {
		databuf.intbuf[0] |= (1L << j);
	      } else {
		databuf.intbuf[0] &= ~(1L << j);
	      }
	      break;
	    }
	    case TYPE_DOUBLE:
	      if (Tcl_GetDouble(interp, argv[2 + j], databuf.doublebuf + j))
		return (TCL_ERROR);
	      break;
	    default: ;
	    }
	  }

	  if (find_callback(i)(interp, databuf.intbuf) != TCL_OK)
	    return gather_runtime_errors(interp, TCL_ERROR);
	  /* fall through to write out the set value immediately again */
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
	case TYPE_BOOL: {
	  if ((*(int *)fields[i].data) & (1L << j))
	    strcpy(buffer, "1");
	  else
	    strcpy(buffer, "0");
	  break;
	}
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
      if (!all) {
	if (writing)
	  return gather_runtime_errors(interp, TCL_OK);
	else
	  return (TCL_OK);
      }
    }
  }
  if (all)
    return TCL_OK;

  Tcl_AppendResult(interp, "unknown md variable \"",
		   argv[1], "\"", (char *) NULL);
  return (TCL_ERROR);
}
