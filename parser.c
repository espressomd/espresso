// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
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
