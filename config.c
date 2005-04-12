// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.

/** \file config.c
 *
 *  contains \ref code_info and version stuff.
 *
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
*/
#include <tcl.h>
#include "utils.h"

int version_callback(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, PACKAGE_NAME ": " PACKAGE_VERSION ", Last Change: " LAST_CHANGE, (char *) NULL);
  return (TCL_OK);
}

/** callback for compilation status. */
int compilation_callback(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "{ Compilation status ", (char *) NULL);
  Tcl_AppendResult(interp, "{ " COMPILE_MODE " } ", (char *) NULL);
  Tcl_AppendResult(interp, "{ MPI " MPI " } ", (char *) NULL);
#ifdef USEFFTW3
  Tcl_AppendResult(interp, "{ FFTW3 } ", (char *) NULL);
#else
  Tcl_AppendResult(interp, "{ FFTW2 } ", (char *) NULL);
#endif
#ifdef TK
  Tcl_AppendResult(interp, "{ TK } ", (char *) NULL);
#endif
#ifdef PARTIAL_PERIODIC
  Tcl_AppendResult(interp, "{ PARTIAL_PERIODIC } ", (char *) NULL);
#endif
#ifdef LJ_WARN_WHEN_CLOSE
  Tcl_AppendResult(interp, "{ LJ_WARN_WHEN_CLOSE } ", (char *) NULL);
#endif
#ifdef ELECTROSTATICS
  Tcl_AppendResult(interp, "{ ELECTROSTATICS } ", (char *) NULL);
#endif
#ifdef ROTATION
  Tcl_AppendResult(interp, "{ ROTATION } ", (char *) NULL);
#endif
#ifdef MASS
  Tcl_AppendResult(interp, "{ MASS } ", (char *) NULL);
#endif
#ifdef EXTERNAL_FORCES
  Tcl_AppendResult(interp, "{ EXTERNAL_FORCES } ", (char *) NULL);
#endif
#ifdef CONSTRAINTS
  Tcl_AppendResult(interp, "{ CONSTRAINTS } ", (char *) NULL);
#endif
#ifdef BOND_CONSTRAINT
  Tcl_AppendResult(interp, "{ BOND_CONSTRAINT } ", (char *) NULL);
#endif
#ifdef EXCLUSIONS
  Tcl_AppendResult(interp, "{ EXCLUSIONS } ", (char *) NULL);
#endif
#ifdef COMFORCE
  Tcl_AppendResult(interp, "{ COMFORCE } ", (char *) NULL);
#endif
#ifdef COMFIXED
  Tcl_AppendResult(interp, "{ COMFIXED } ", (char *) NULL);
#endif
#ifdef TABULATED
  Tcl_AppendResult(interp, "{ TABULATED } ", (char *) NULL);
#endif
#ifdef LENNARD_JONES
  Tcl_AppendResult(interp, "{ LENNARD_JONES } ", (char *) NULL);
#endif
#ifdef MORSE
  Tcl_AppendResult(interp, "{ MORSE } ", (char *) NULL);
#endif
#ifdef BUCKINGHAM
  Tcl_AppendResult(interp, "{ BUCKINGHAM } ", (char *) NULL);
#endif
#ifdef SOFT_SPHERE
  Tcl_AppendResult(interp, "{ SOFT_SPHERE } ", (char *) NULL);
#endif
#ifdef LJCOS
  Tcl_AppendResult(interp, "{ LJCOS } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_HARMONIC
  Tcl_AppendResult(interp, "{ BOND_ANGLE_HARMONIC } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_COSINE
  Tcl_AppendResult(interp, "{ BOND_ANGLE_COSINE } ", (char *) NULL);
#endif
#ifdef BOND_ANGLE_COSSQUARE
  Tcl_AppendResult(interp, "{ BOND_ANGLE_COSSQUARE } ", (char *) NULL);
#endif
#ifdef NEMD
  Tcl_AppendResult(interp, "{ NEMD } ", (char *) NULL);
#endif
#ifdef NPT
  Tcl_AppendResult(interp, "{ NPT } ", (char *) NULL);
#endif
#ifdef DPD
  Tcl_AppendResult(interp, "{ DPD } ", (char *) NULL);
#endif
#ifdef LB
  Tcl_AppendResult(interp, "{ LB } ", (char *) NULL);
#endif
  Tcl_AppendResult(interp, "}", (char *) NULL);
  return (TCL_OK);
}
