#ifndef EXTERNAL_POTENTIAL_TCL_H
#define EXTERNAL_POTENTIAL_TCL_H


#include <tcl.h>
#include "parser.hpp"
#include "external_potential.hpp"

int tclcommand_external_potential(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv);

#endif
