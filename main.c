#include <tcl.h>
#include "initialize.h"
#include <mpi.h>
#include "global.h"
#include "slave.h"

/* initialize MPI and determine nprocs/node */
void init_mpi(int *argc, char ***argv)
{
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &node);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
}

int appinit(Tcl_Interp *interp)
{
  if (Tcl_Init(interp) == TCL_ERROR)
    return (TCL_ERROR);
#ifdef USE_TK
  if (Tk_Init(interp) == TCL_ERROR)
    return (TCL_ERROR);
#endif
  if (initialize(interp) == TCL_ERROR)
    return (TCL_ERROR);
  return (TCL_OK);
}

int main(int argc, char **argv)
{
  init_mpi(&argc, &argv);

  if (node == 0) {
    /* master node */
    Tcl_Main(argc, argv, appinit);
    return 0;
  }

  /* slave node */
  process_loop();

  return 0;
}
