// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
/** \file main.c
    Main file of Espresso. Initialization of tcl interpreter and exit handling.

    DO NOT CHANGE!!!

*/
/* first, since we need the TK define */
#include "utils.h"
#include <tcl.h>
#ifdef TK
#include <tk.h>
#endif
#include <stdlib.h>
#include <signal.h>
#include "initialize.h"
#include "global.h"
#include "communication.h"

void errexit()
{
#ifdef FORCE_CORE
  core();
#endif
  exit(1);
}

void sigint_handler(int sig) 
{
  /* without this exit handler the nodes might exit asynchronously
   * without calling MPI_Finalize, which may cause MPI to hang
   * (e.g. this handler makes CTRL-c work properly with poe)
   *
   * NOTE: mpirun installs its own handler for SIGINT/SIGTERM
   *       and takes care of proper cleanup and exit */

  char *errtxt;

  static int numcalls = 0;
  if (numcalls++ > 0) exit(sig); // catch sig only once

  /* we use runtime_error to indicate that sig was called;
   * upon next call of mpi_gather_runtime_errors all nodes 
   * will clean up and exit. */
  errtxt = runtime_error(64);
  ERROR_SPRINTF(errtxt, "{000 caught signal %d} ",sig);
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
#ifdef EFENCE
  extern int EF_ALLOW_MALLOC_0;
  EF_ALLOW_MALLOC_0 = 1;
#endif

  mpi_init(&argc, &argv);

  /* register handler for SIGINT */
  signal(SIGINT, sigint_handler);

  if (this_node == 0) {
    /* master node */
#ifdef FORCE_CORE
    /* core should be the last exit handler (process dies) */
    atexit(core);
#endif
    atexit(mpi_stop);
#ifdef TK
    Tk_Main(argc, argv, appinit);
#else
    Tcl_Main(argc, argv, appinit);
#endif
  } 
  else {
    /* slave node */
    on_program_start(0);
    mpi_loop();
  }

  return 0;
}
