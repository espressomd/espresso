#ifndef PARSER_H
#define PARSER_H
/** \file parser.h
    This file contains macros for parsing the parameters to the
    'inter' command.

    <b>Responsible:</b>
    <a href="mailto:schneidc@mpip-mainz.mpg.de">Christoph</a>
 */

#include <string.h>

#define ARG_IS_S(no, str) !strncasecmp(argv[(no)], (str), strlen(argv[(no)]))
#define ARG0_IS_S(str) ARG_IS_S(0, (str))
#define ARG1_IS_S(str) ARG_IS_S(1, (str))

#define ARG_IS_I(no, dest) (!(Tcl_GetInt(interp, argv[(no)], &(dest)) == TCL_ERROR))
#define ARG0_IS_I(dest) ARG_IS_I(0, (dest))
#define ARG1_IS_I(dest) ARG_IS_I(1, (dest))

#define ARG_IS_D(no, dest) (!(Tcl_GetDouble(interp, argv[(no)], &(dest)) == TCL_ERROR))
#define ARG0_IS_D(dest) ARG_IS_D(0, (dest))
#define ARG1_IS_D(dest) ARG_IS_D(1, (dest))

#define CHECK_VALUE(func, errmsg) \
        if (func == TCL_ERROR) { \
           Tcl_AppendResult(interp, errmsg, (char *)NULL); \
           return TCL_ERROR; \
        } else \
	  return TCL_OK;
#endif
