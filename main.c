/** \file main.c
    Main file of tcl_md. Initialization of tcl interpreter and exit handling.

    DO NOT CHANGE!!!

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

*/
#include <tcl.h>
#ifdef TK
#include <tk.h>
#endif
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

#ifdef TK
  if (Tk_Init(interp) == TCL_ERROR)
    return (TCL_ERROR);
#endif

  if (on_program_start(interp) == TCL_ERROR)
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
#ifdef TK
    Tk_Main(argc, argv, appinit);
#else
    Tcl_Main(argc, argv, appinit);
#endif
    return 0;
  }
  else
    on_program_start(0);

  /* slave node */
  mpi_loop();

  return 0;
}
