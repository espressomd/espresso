#ifndef PARSER_H
#define PARSER_H
/** \file parser.h
    This file contains macros for parsing the parameters to the
    'inter' command.

    <b>Responsible:</b>
    <a href="mailto:schneidc@mpip-mainz.mpg.de">Christoph</a>
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

#define CHECK_VALUE(func, errmsg) \
        if (func == TCL_ERROR) { \
           Tcl_AppendResult(interp, errmsg, (char *)NULL); \
           return TCL_ERROR; \
        } else \
	  return TCL_OK;
#endif
