// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
/** \file adresso.c
    This is the place for adaptive resolution scheme
    Implementation of adresso.h
*/

#include "adresso.h"
#include "communication.h"
#include "parser.h"
#include "cells.h"

/** \name Privat Functions */
/************************************************************/
/*@{*/

/*@}*/

double adress_vars[7]       = {0, 0, 0, 0, 0, 0, 0};

int adress_tcl(ClientData data, Tcl_Interp *interp, int argc, char **argv){
   int err = TCL_OK;
   Tcl_ResetResult(interp);
   Tcl_AppendResult(interp, "Adress is not compiled in (change config.h).", (char *)NULL);
   err = (TCL_ERROR);
   return mpi_gather_runtime_errors(interp, err);
}

