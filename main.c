/** \file main.c
    Main file of tcl_md. Initialization of tcl interpreter and exit handling.
*/
#include <tcl.h>
#include <stdlib.h>
#include "initialize.h"
#include "global.h"
#include "communication.h"
#include "debug.h"

void errexit()
{
#ifdef FORCE_CORE
  core();
#endif
  exit(1);
}

void exitHandler(ClientData data)
{
  mpi_stop();
}

int appinit(Tcl_Interp *interp)
{
  if (Tcl_Init(interp) == TCL_ERROR)
    return (TCL_ERROR);
  Tcl_CreateExitHandler(exitHandler, 0);

  if (initialize(interp) == TCL_ERROR)
    return (TCL_ERROR);
  return (TCL_OK);
}

int main(int argc, char **argv)
{
  mpi_init(&argc, &argv);

  if (this_node == 0) {
    /* master node */
    atexit(mpi_stop);
#ifdef FORCE_CORE
    atexit(core);
#endif
    Tcl_Main(argc, argv, appinit);
    return 0;
  }
  else
    initialize(0);

  /* slave node */
  mpi_loop();

  return 0;
}
