#include "initialize.h"
#include "field.h"
#include "chain.h"

int initialize(Tcl_Interp *interp)
{
  /*
    call the initialization of the modules here
    i. e. adding the commands etc.
  */

  /* field.c */
  field_init(interp);

  /* chain.c */
  chain_init(interp);

  return TCL_OK;
}
