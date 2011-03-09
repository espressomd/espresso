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
#ifndef PARSER_H
#define PARSER_H
/** \file parser.h
    This file contains macros for parsing the parameters to the
    'inter' command.
 */

#include <string.h>
#include <tcl.h>
#include "utils.h"
/** parse an integer list
    @param interp for conversion of backslash stuff
    @param list the string containing the list
    @param il where to store the results
*/
int parse_int_list(Tcl_Interp *interp, char *list, IntList *il);

/** parse an double list
    @param interp for conversion of backslash stuff
    @param list the string containing the list
    @param dl where to store the results
*/
int parse_double_list(Tcl_Interp *interp, char *list, DoubleList *dl);


#define ARG_IS_S(no, str) !strncasecmp(argv[(no)], (str), strlen(argv[(no)]))
#define ARG0_IS_S(str) ARG_IS_S(0, (str))
#define ARG1_IS_S(str) ARG_IS_S(1, (str))

#define ARG_IS_I(no, dest) (!(Tcl_GetInt(interp, argv[(no)], &(dest)) == TCL_ERROR))
#define ARG0_IS_I(dest) ARG_IS_I(0, (dest))
#define ARG1_IS_I(dest) ARG_IS_I(1, (dest))

#define ARG_IS_D(no, dest) (!(Tcl_GetDouble(interp, argv[(no)], &(dest)) == TCL_ERROR))
#define ARG0_IS_D(dest) ARG_IS_D(0, (dest))
#define ARG1_IS_D(dest) ARG_IS_D(1, (dest))

#define ARG_IS_INTLIST(no, il) (parse_int_list(interp, argv[(no)], &(il)))
#define ARG0_IS_INTLIST(il) ARG_IS_INTLIST(0, (il))
#define ARG1_IS_INTLIST(il) ARG_IS_INTLIST(1, (il))

#define ARG_IS_DOUBLELIST(no, il) (parse_double_list(interp, argv[(no)], &(il)))
#define ARG0_IS_DOUBLELIST(il) ARG_IS_DOUBLELIST(0, (il))
#define ARG1_IS_DOUBLELIST(il) ARG_IS_DOUBLELIST(1, (il))

#define CHECK_VALUE(func, errmsg) \
        if (func == TCL_ERROR) { \
           Tcl_AppendResult(interp, errmsg, (char *)NULL); \
           return TCL_ERROR; \
        } \
	return TCL_OK;
#endif
