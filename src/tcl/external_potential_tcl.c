

#include "external_potential_tcl.h"
int tclcommand_external_potential(ClientData _data, Tcl_Interp *interp,
	       int argc, char **argv) {
  ExternalPotential* e;
  int error = generate_external_potential(&e);
  if (error == ES_ERROR)
    return TCL_ERROR;
  else
    return TCL_OK;

}
