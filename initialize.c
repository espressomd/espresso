#include "initialize.h"
#include "tcl_datafield.h"

int initialize(Tcl_Interp *interp)
{
  /*
    call the initialization of the modules here
    i. e. adding the commands etc.
  */

  tcl_datafield_init(interp);

  return TCL_OK;
}
