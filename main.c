#include <tcl.h>
#include "initialize.h"
#include "global.h"
#include "communication.h"

void exitHandler(ClientData data)
{
  stop_mpi();
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
  init_mpi(&argc, &argv);

  if (this_node == 0) {
    /* master node */
    init_data();
    Tcl_Main(argc, argv, appinit);
    return 0;
  }

  /* slave node */
  init_data();
  mpi_loop();

  return 0;
}
