// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
/** \file main.c
    Main file of Espresso. Initialization of tcl interpreter and exit handling.

    DO NOT CHANGE!!!

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

*/
/* first, since we need the TK define */
#include "utils.h"
#include <tcl.h>
#ifdef TK
#include <tk.h>
#endif
#include <stdlib.h>
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
