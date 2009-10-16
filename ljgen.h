// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2007; all rights reserved unless otherwise stated.
#ifndef LJGEN_H
#define LJGEN_H

/** \file ljgen.h
 *  Routines to calculate the generalized lennard jones energy and/or  force 
 *  for a particle pair. "Generalized" here means that the LJ energy is of the
 *  form
 *
 *  eps * [ b1 * (sigma/(r-r_offset))^a1 - b2 * (sigma/(r-r_offset))^a2 + shift]
 *
 *  \ref forces.c
*/

/* These headers are needed to define types used in this header, hence they
 * are included here.  */
#include "particle_data.h"
#include "interaction_data.h"

int printljgenIAToResult(Tcl_Interp *interp, int i, int j);

int ljgen_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv);


/** Calculate lennard Jones force between particle p1 and p2 */
void add_ljgen_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3]);

/** calculate Lennard jones energy between particle p1 and p2. */
double ljgen_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist);

/** calculate lj_capradius from lj_force_cap */
void calc_ljgen_cap_radii(double force_cap);

/* LJGEN_H */
#endif 
