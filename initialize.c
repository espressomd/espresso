#include "initialize.h"
#include "tcl_datafield.h"
#include "integrate.h"

int initialize(Tcl_Interp *interp)
{
  /*
    call the initialization of the modules here
    i. e. adding the commands etc.
  */

  tcl_datafield_init(interp);

  tcl_integrator_init(interp);

  return TCL_OK;
}
