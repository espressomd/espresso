#include <tcl.h>
#include <stdlib.h>
#include "initialize.h"
#include "global.h"
#include "communication.h"

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

    init_data();
    Tcl_Main(argc, argv, appinit);
    return 0;
  }

  /* slave node */
  init_data();
  mpi_loop();

  return 0;
}
