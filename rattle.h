/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef RATTLE_H
#define RATTLE_H

/** \file rattle.h    RATTLE Algorithm (Rattle: A "Velocity" Version of the Shake
 *                    Algorithm for Molecular Dynamics Calculations, H.C Andersen,
 *                    J Comp Phys, 52, 24-34, 1983)
 *
 *  For more information see \ref rattle.c "rattle.c".
*/
#include "parser.h"
#include "global.h"
#include "particle_data.h"
#include "integrate.h"

/** number of rigid bonds */
extern int n_rigidbonds;

#ifdef BOND_CONSTRAINT

/** Transfers the current particle positions from r.p[3] to r.p_pold[3]
    of the \ref Particle structure. Invoked from \ref correct_pos_shake() */
void save_old_pos();

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

  if (time_step < 0. ) {
    Tcl_AppendResult(interp, "set time_step before declaring rigid_bond" , (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(rigid_bond_set_params(bond_type, d, p_tol, v_tol), "bond type must be nonnegative");
}
#endif
#endif
