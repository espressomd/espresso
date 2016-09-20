/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _TCL_PARSER_H
#define _TCL_PARSER_H
/** \file parser.hpp
    This file contains macros for parsing the parameters to the
    'inter' command.
 */

#include "utils.hpp"
#include <cstring>
#include <tcl.h>

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

/** gather all error messages from all nodes and set the interpreter result
    to these error messages. This should be called only on the master node.

    The errors are append to the result, if ret_state == TCL_ERROR,
    otherwise the result is overwritten, in case an error occurs.
    Therefore you should end any Tcl command handler by return
    gather_runtime_errors(<return_value>).  This code uses
    asynchronous communication.

    @param ret_state return value of the procedure
    @param interp where to put the errors
    @return new return value after the background errors, if any, have been handled
*/
int gather_runtime_errors(Tcl_Interp *interp, int ret_state);


#define ARG_IS_S(no, str) !strncasecmp(argv[(no)], (str), strlen(argv[(no)]))
#define ARG0_IS_S(str) ARG_IS_S(0, (str))
#define ARG1_IS_S(str) ARG_IS_S(1, (str))
#define ARG_IS_S_EXACT(no, str) !strcmp(argv[(no)], (str))
#define ARG0_IS_S_EXACT(str) ARG_IS_S_EXACT(0, (str))
#define ARG1_IS_S_EXACT(str) ARG_IS_S_EXACT(1, (str))

#define ARG_IS_I(no, dest) (!(Tcl_GetInt(interp, argv[(no)], &(dest)) == TCL_ERROR))
#define ARG0_IS_I(dest) ARG_IS_I(0, (dest))
#define ARG1_IS_I(dest) ARG_IS_I(1, (dest))

#define ARG_IS_D(no, dest) ((argc > 0) && (!(Tcl_GetDouble(interp, argv[(no)], &(dest)) == TCL_ERROR)))
#define ARG0_IS_D(dest) ARG_IS_D(0, (dest))
#define ARG1_IS_D(dest) ARG_IS_D(1, (dest))

#define ARG_IS_INTLIST(no, il) (parse_int_list(interp, argv[(no)], &(il)))
#define ARG0_IS_INTLIST(il) ARG_IS_INTLIST(0, (il))
#define ARG1_IS_INTLIST(il) ARG_IS_INTLIST(1, (il))

#define ARG_IS_DOUBLELIST(no, il) (parse_double_list(interp, argv[(no)], &(il)))
#define ARG0_IS_DOUBLELIST(il) ARG_IS_DOUBLELIST(0, (il))
#define ARG1_IS_DOUBLELIST(il) ARG_IS_DOUBLELIST(1, (il))

#define CHECK_VALUE(func, errmsg) \
        if (func == ES_ERROR) { \
           Tcl_AppendResult(interp, errmsg, (char *)NULL); \
           return TCL_ERROR; \
        } \
	return TCL_OK;
#endif
