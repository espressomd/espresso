// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

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

#ifdef BOND_CONSTRAINT
#ifdef ROTATION
#error BOND_CONSTRAINT and ROTATION currently do not work together
#endif
#endif

/* errors for all modules that require fftw if not present */
#ifndef FFTW

#ifdef MODES
#error MODES requires the fftw
#endif

#ifdef LB
#error LB requires the fftw
#endif

#endif

int version_callback(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, PACKAGE_NAME ": " PACKAGE_VERSION ", Last Change: " LAST_CHANGE, (char *) NULL);
  return (TCL_OK);
}

/** callback for compilation status. */
int compilation_callback(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "{ Compilation status ", (char *) NULL);
#ifdef DEBUG
  Tcl_AppendResult(interp, "{ " DEBUG " } ", (char *) NULL);
#endif
#ifdef PROFILING
  Tcl_AppendResult(interp, "{ " PROFILING " } ", (char *) NULL);
#endif
  Tcl_AppendResult(interp, "{ MPI " MPI " } ", (char *) NULL);
#if FFTW == 2
  Tcl_AppendResult(interp, "{ FFTW2 } ", (char *) NULL);
#elif FFTW == 3
  Tcl_AppendResult(interp, "{ FFTW3 } ", (char *) NULL);
#endif
#ifdef TK
  Tcl_AppendResult(interp, "{ TK } ", (char *) NULL);
#endif
#ifdef MODES
  Tcl_AppendResult(interp, "{ MODES } ", (char *) NULL);
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
#ifdef DIPOLES
  Tcl_AppendResult(interp, "{ DIPOLES } ", (char *) NULL);
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
#ifdef BOND_VIRTUAL
  Tcl_AppendResult(interp, "{ BOND_VIRTUAL } ", (char *) NULL);
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
#ifdef LENNARD_JONES_GENERIC
  Tcl_AppendResult(interp, "{ LENNARD_JONES_GENERIC } ", (char *) NULL);
#endif
#ifdef LJ_ANGLE
  Tcl_AppendResult(interp, "{ LJ_ANGLE } ", (char *) NULL);
#endif
#ifdef SMOOTH_STEP
  Tcl_AppendResult(interp, "{ SMOOTH_STEP } ", (char *) NULL);
#endif
#ifdef BMHTF_NACL
  Tcl_AppendResult(interp, "{ BMHTF_NACL } ", (char *) NULL);
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
#ifdef LJCOS2
  Tcl_AppendResult(interp, "{ LJCOS2 } ", (char *) NULL);
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
#ifdef TRANS_DPD
  Tcl_AppendResult(interp, "{ TRANS_DPD } ", (char *) NULL);
#endif
#ifdef DPD
  Tcl_AppendResult(interp, "{ DPD } ", (char *) NULL);
#endif
#ifdef DPD_MASS
  Tcl_AppendResult(interp, "{ DPD_MASS } ", (char *) NULL);
#endif
#ifdef LB
  Tcl_AppendResult(interp, "{ LB } ", (char *) NULL);
#endif
#ifdef INTER_DPD
  Tcl_AppendResult(interp, "{ INTER_DPD } ", (char *) NULL);
#endif
#ifdef INTER_RF
  Tcl_AppendResult(interp, "{ INTER_RF } ", (char *) NULL);
#endif
#ifdef NO_INTRA_NB
  Tcl_AppendResult(interp, "{ NO_INTRA_NB } ", (char *) NULL);
#endif
  Tcl_AppendResult(interp, "}", (char *) NULL);
  return (TCL_OK);
}
