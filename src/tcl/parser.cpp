/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file parser.cpp
    Implementation of \ref parser.hpp "parser.h". \ref parse_int_list is too long for inlining.
 */

#include "utils.hpp"
#include "parser.hpp"
#include "communication.hpp"

int parse_int_list(Tcl_Interp *interp, char *list, IntList *il)
{
  int i, tmp_argc, res = 1;
  char  **tmp_argv;
  Tcl_SplitList(interp, list, &tmp_argc, &tmp_argv);
  realloc_intlist(il, il->n = tmp_argc);
  for(i = 0 ; i < tmp_argc; i++) if (Tcl_GetInt(interp, tmp_argv[i], &(il->e[i])) == TCL_ERROR) { res = 0; break; }
  Tcl_Free((char *)tmp_argv);
  return res;
}

int parse_double_list(Tcl_Interp *interp, char *list, DoubleList *dl)
{
  int i, tmp_argc, res = 1;
  char  **tmp_argv;
  Tcl_SplitList(interp, list, &tmp_argc, &tmp_argv);
  realloc_doublelist(dl, dl->n = tmp_argc);
  for(i = 0 ; i < tmp_argc; i++) if (Tcl_GetDouble(interp, tmp_argv[i], &(dl->e[i])) == TCL_ERROR) { res = 0; break; }
  Tcl_Free((char *)tmp_argv);
  return res;
}

int gather_runtime_errors(Tcl_Interp *interp, int error_code)
{
  char **errors =(char **)malloc(n_nodes*sizeof(char *));

  if (mpi_gather_runtime_errors(errors) == ES_OK) {
    free(errors);
    return error_code;
  }

  /* reset any results of the previous command, since we got an error
     during evaluation, they are at best bogus. But any normal error
     messages should be kept, they might help to track down the
     problem.
  */
  if (error_code != TCL_ERROR)
    Tcl_ResetResult(interp);
  else
    Tcl_AppendResult(interp, " ", (char *) NULL);

  Tcl_AppendResult(interp, "background_errors ", (char *) NULL);

  for (int node = 0; node < n_nodes; node++) {
    if (errors[node] != NULL) {
      char nr_buf[TCL_INTEGER_SPACE + 3];
      sprintf(nr_buf, "%d ", node);
      
      /* check whether it's the same message as from the master, then
	 just consent */
      if (node > 0 && errors[0] && strcmp(errors[node], errors[0]) == 0)
	Tcl_AppendResult(interp, nr_buf, "<consent> ", (char *) NULL);
      else
	Tcl_AppendResult(interp, nr_buf, errors[node], (char *) NULL);

      /* we still need the master nodes error for squashing, see
	 above */
      if (node > 0)
	free(errors[node]);
    }
  }

  free(errors[0]);
  free(errors);

  return TCL_ERROR;
}
