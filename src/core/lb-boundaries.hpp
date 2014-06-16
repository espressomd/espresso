/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file lb-boundaries.hpp
 *
 * Boundary conditions for Lattice Boltzmann fluid dynamics.
 * Header file for \ref lb-boundaries.cpp.
 *
 * In the current version only simple bounce back walls are implemented. Thus
 * after the streaming step, in all wall nodes all populations are bounced
 * back from where they came from. Ulf Schiller spent a lot of time
 * working on more powerful alternatives, they are to be found in the
 * lb_testing branch of espresso until the end of 2010. Now we stripped
 * down the code to a minimum, as most of it was not sufficiently understandable.
 *
 * Anyone who wants to revive these, please look into the git.
 *
 */

#ifndef LB_BOUNDARIES_H
#define LB_BOUNDARIES_H

#include "utils.hpp"
#include "constraint.hpp"

#if defined (LB_BOUNDARIES) || defined (LB_BOUNDARIES_GPU)

/** wall constraint applied */
#define LB_BOUNDARY_WAL 1
/** spherical constraint applied */
#define LB_BOUNDARY_SPH 2
/** (finite) cylinder shaped constraint applied */
#define LB_BOUNDARY_CYL 3
/** a pore geometry */
#define LB_BOUNDARY_POR 4
/** rhomboid shaped constraint applied */
#define LB_BOUNDARY_RHOMBOID 5
/** stomatocyte shaped constraint applied */
#define LB_BOUNDARY_STOMATOCYTE 6
/** box shaped constraint (for electrokinetics reaction code) */
#define LB_BOUNDARY_BOX 7
/** hollow cone shaped constraint (for electrokinetics reaction code) */
#define LB_BOUNDARY_HOLLOW_CONE 8

// If we have several possible types of boundary treatment
#define LB_BOUNDARY_BOUNCE_BACK 1

/** Structure to specify a boundary. */
typedef struct {
  /** type of the boundary. */
  int type;
  double slip_pref;

  union {
    Constraint_wall wal;
    Constraint_sphere sph;
    Constraint_cylinder cyl;
    Constraint_rhomboid rhomboid;
    Constraint_pore pore;
    Constraint_stomatocyte stomatocyte;
    Constraint_box box;
    Constraint_hollow_cone hollow_cone;
  } c;
  double force[3];
  double velocity[3];
#ifdef EK_BOUNDARIES
  float charge_density;
#endif
} LB_Boundary;

extern int n_lb_boundaries;
extern LB_Boundary *lb_boundaries;

/*@}*/

/** Initializes the constrains in the system. 
 *  This function determines the lattice sited which belong to boundaries
 *  and marks them with a corresponding flag. 
 */
void lb_init_boundaries();
void lbboundary_mindist_position(double pos[3], double* mindist, double distvec[3], int* no); 

int lbboundary_get_force(int no, double* f); 

#endif // (LB_BOUNDARIES) || (LB_BOUNDARIES_GPU)

#ifdef LB_BOUNDARIES
/** Bounce back boundary conditions.
 * The populations that have propagated into a boundary node
 * are bounced back to the node they came from. This results
 * in no slip boundary conditions.
 *
 * [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 */
void lb_bounce_back();
#endif /* LB_BOUNDARIES */

#endif /* LB_BOUNDARIES_H */
