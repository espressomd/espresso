#include <tcl.h>
#include <stdlib.h>
#include "initialize.h"
#include "global.h"
#include "communication.h"
#include "debug.h"

#ifdef FORCE_CORE
int regular_exit = 0;
static int core_done = 0;

void core()
{
  if (!core_done && !regular_exit) {
    core_done = 1;
    fprintf(stderr, "forcing core dump on irregular exit\n");
    *(int *)0 = 0;
  }
}
#endif

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

  /* slave node */
  mpi_loop();

  return 0;
}
