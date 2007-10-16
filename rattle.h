// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef RATTLE_H
#define RATTLE_H

/** \file rattle.h    RATTLE Algorithm (Rattle: A "Velocity" Version of the Shake
 *                    Algorithm for Molecular Dynamics Calculations, H.C Andersen,
 *                    J Comp Phys, 52, 24-34, 1983)
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arijitmaitra@uni-muenster.de">Arijit</a>
 *
 *  For more information see \ref rattle.c "rattle.c".
*/
#include "parser.h"
#include "global.h"
#include "particle_data.h"

/** number of rigid bonds */
extern int n_rigidbonds;

#ifdef BOND_CONSTRAINT

/** Transfers the current particle positions from r.p[3] to r.p_pold[3]
    of the \ref Particle structure. Invoked from \ref correct_pos_shake() */
void save_old_pos();
void init_f_shake();
/** Propagate velocity and position while using SHAKE algorithm for bond constraint.*/
void correct_pos_shake();

/** Correction of current velocities using RATTLE algorithm*/
void correct_vel_shake();

/** set the parameter for a rigid, aka RATTLE bond */
int rigid_bond_set_params(int bond_type, double d, double p_tol, double v_tol);

/// parse parameters for the rigid bonds
MDINLINE int inter_parse_rigid_bonds(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double d, p_tol, v_tol;
    
  if (argc != 4) {
    Tcl_AppendResult(interp, "rigid bond needs 3 parameters: "
		     "<constrained_bond_distance> <Positional_tolerance> <Velocity_tolerance>", (char *) NULL);
    return TCL_ERROR;
  }

  if ((! ARG_IS_D(1, d)) || (! ARG_IS_D(2, p_tol)) || (! ARG_IS_D(3, v_tol)) ) {
    Tcl_AppendResult(interp, "rigid bond needs 3 DOUBLE parameters: "
		     "<constrained_bond_distance> <Positional_tolerance> <Velocity_tolerance>", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(rigid_bond_set_params(bond_type, d, p_tol, v_tol), "bond type must be nonnegative");
}
#endif
#endif
