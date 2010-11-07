/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
/** \file parser.c
    Implementation of \ref parser.h "parser.h". \ref parse_int_list is too long for inlining.
 */

#include "utils.h"
#include "parser.h"

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
