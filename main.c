#include <tcl.h>
#include <stdlib.h>
#include "initialize.h"
#include "global.h"
#include "communication.h"

#ifdef FORCE_CORE
void core()
{
  *(int *)0 = 0;
}
#endif

void exitHandler(ClientData data)
{
  mpi_stop();
}

#ifdef DEFINE_TCL_RAND
int tcl_rand(ClientData clientData, Tcl_Interp *interp,
	     Tcl_Value *args, Tcl_Value *resultPtr)
{
  resultPtr->type = TCL_DOUBLE;
  resultPtr->doubleValue = drand48();
  return TCL_OK;
}

int tcl_srand(ClientData clientData, Tcl_Interp *interp,
	      Tcl_Value *args, Tcl_Value *resultPtr)
{
  int i;
  srand48(args[0].intValue);
  return tcl_rand(0, interp, NULL, resultPtr);
}

#endif

int appinit(Tcl_Interp *interp)
{
#ifdef DEFINE_TCL_RAND
  Tcl_ValueType vt[2] = { TCL_INT };
#endif

  if (Tcl_Init(interp) == TCL_ERROR)
    return (TCL_ERROR);
  Tcl_CreateExitHandler(exitHandler, 0);

#ifdef DEFINE_TCL_RAND
  Tcl_CreateMathFunc(interp, "rand", 0, NULL, tcl_rand, 0);
  Tcl_CreateMathFunc(interp, "srand", 1, vt, tcl_srand, 0);
#endif

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
