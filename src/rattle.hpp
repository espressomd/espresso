/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

/** \file rattle.hpp    RATTLE Algorithm (Rattle: A "Velocity" Version of the Shake
 *                    Algorithm for Molecular Dynamics Calculations, H.C Andersen,
 *                    J Comp Phys, 52, 24-34, 1983)
 *
 *  For more information see \ref rattle.cpp "rattle.c".
*/
#include "global.hpp"
#include "particle_data.hpp"
#include "integrate.hpp"

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

#endif
#endif
