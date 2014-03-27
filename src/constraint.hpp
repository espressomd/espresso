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
#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

/** \file constraint.hpp
 *  Routines for handling of constraints.
 *  Only active if the feature CONSTRAINTS is activated.
 *  see also \ref interaction_data.hpp
 */

#include "particle_data.hpp"
#include "interaction_data.hpp"

#ifdef CONSTRAINTS

/** number of constraints. */
extern int n_constraints;
/** field containing constraints. */
extern Constraint *constraints;

extern int reflection_happened;

/** Exported functions
 */

Constraint *generate_constraint();

void calculate_wall_dist(Particle *p1, double ppos[3], 
			 Particle *c_p, Constraint_wall *c, 
			 double *dist, double *vec);

void calculate_sphere_dist(Particle *p1, double ppos[3], 
			   Particle *c_p, Constraint_sphere *c, 
			   double *dist, double *vec);

void calculate_maze_dist(Particle *p1, double ppos[3], 
			 Particle *c_p, Constraint_maze *c, 
			 double *dist, double *vec);

void calculate_cylinder_dist(Particle *p1, double ppos[3], 
			     Particle *c_p, Constraint_cylinder *c, 
			     double *dist, double *vec);

void calculate_rhomboid_dist(Particle *p1, double ppos[3], 
			     Particle *c_p, Constraint_rhomboid *c, 
			     double *dist, double *vec);

void calculate_pore_dist(Particle *p1, double ppos[3], 
			 Particle *c_p, Constraint_pore *c, 
			 double *dist, double *vec);

void calculate_slitpore_dist(Particle *p1, double ppos[3], 
			 Particle *c_p, Constraint_slitpore *c, 
			 double *dist, double *vec);

void calculate_plane_dist(Particle *p1, double ppos[3], 
			  Particle *c_p, Constraint_plane *c, 
			  double *dist, double *vec);

void calculate_stomatocyte_dist( Particle *p1, double ppos [3], 
        Particle *c_p, Constraint_stomatocyte *cons, 
        double *dist, double *vec );

void calculate_hollow_cone_dist( Particle *p1, double ppos [3], 
        Particle *c_p, Constraint_hollow_cone *cons, 
        double *dist, double *vec );

void add_rod_force(Particle *p1, double ppos[3], 
		   Particle *c_p, Constraint_rod *c);

double rod_energy(Particle *p1, double ppos[3], 
		  Particle *c_p, Constraint_rod *c);

void add_plate_force(Particle *p1, double ppos[3], 
		     Particle *c_p, Constraint_plate *c);

double plate_energy(Particle *p1, double ppos[3], 
		    Particle *c_p, Constraint_plate *c);

void add_ext_magn_field_force(Particle *p1, 
			      Constraint_ext_magn_field *c);

double ext_magn_field_energy(Particle *p1, 
			     Constraint_ext_magn_field *c);

void reflect_particle(Particle *p1, double *distance_vec, 
		      int reflecting);

void add_constraints_forces(Particle *p1);

double add_constraints_energy(Particle *p1);

void init_constraint_forces();
#endif

#endif
