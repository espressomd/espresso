// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file forces.c Force calculation.
 *
 *  For more information see \ref forces.h "forces.h".
*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "debug.h"
#include "thermostat.h"
#include "communication.h"
#include "ghosts.h" 
#include "verlet.h"
#include "utils.h"
#include "grid.h"
#include "cells.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "rotation.h"
#include "forces.h"

#if 0

/************************************************************/

/** nonbonded and bonded force calculation using the verlet list */
void calculate_verlet_ia()
{
  Cell *cell;
  Particle *p, **pairs;
  Particle *p1, *p2;
  int k, i,j, c, np;
  double d[3], dist2, dist;
  IA_parameters *ia_params;

  /* force calculation loop. */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->pList.part;
    np = cell->pList.n;

    /* calculate bonded interactions */
    calc_bonded_forces(p, np);

    /* calculate non bonded interactions (loop verlet lists of neighbors) */
    for (k = 0; k < cell->n_neighbors; k++) {
      pairs = cell->nList[k].vList.pair;  /* verlet list */
      np    = cell->nList[k].vList.n;     /* length of verlet list */

      /* verlet list loop */
      for(i=0; i<2*np; i+=2) {
	p1 = pairs[i];                    /* pointer to particle 1 */
	p2 = pairs[i+1];                  /* pointer to particle 2 */

	ia_params = get_ia_param(p1->p.type,p2->p.type);
	/* distance calculation */
	for(j=0; j<3; j++) d[j] = p1->r.p[j] - p2->r.p[j];
	dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
	dist  = sqrt(dist2);

	add_non_bonded_pair_force(p1, p2, ia_params, d, dist, dist2);
      } 
    }
  }
}

void force_calc()
{
  /* preparation */
  init_forces();

  calculate_verlet_ia();

  calc_long_range_forces();
}

/************************************************************/

void calc_bonded_forces(Particle *p, int np)
{
  int j;

  /* calculate bonded interactions (loop local particles) */
  for(j = 0; j < np; j++) {
    add_bonded_pair_force(&p[j]);
#ifdef CONSTRAINTS
    add_constraints_forces(&p[j]);
#endif
  }
}

/************************************************************/

void calc_long_range_forces()
{
#ifdef ELECTROSTATICS  
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
  case COULOMB_P3M:
    P3M_calc_kspace_forces(1,0);
    break;
  case COULOMB_MMM1D:
    MMM1D_calc_forces();
    break;
  }
#endif
}

/************************************************************/

/** initialize the forces for a real particle */
MDINLINE void init_local_particle_force(Particle *part)
{
  friction_thermo(part);
#ifdef EXTERNAL_FORCES   
  if(part->l.ext_flag == PARTICLE_EXT_FORCE) {
    part->f.f[0] += part->l.ext_force[0];
    part->f.f[1] += part->l.ext_force[1];
    part->f.f[2] += part->l.ext_force[2];
  }
#endif

#ifdef ROTATION
  {
    double scale;
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;
    
    /* and rescale quaternion, so it is exactly of unit length */	
    scale = sqrt( SQR(part->r.quat[0]) + SQR(part->r.quat[1]) +
		  SQR(part->r.quat[2]) + SQR(part->r.quat[3]));
    part->r.quat[0]/= scale;
    part->r.quat[1]/= scale;
    part->r.quat[2]/= scale;
    part->r.quat[3]/= scale;
  }
#endif
}

/** initialize the forces for a ghost particle */
MDINLINE void init_ghost_force(Particle *part)
{
  part->f.f[0] = 0;
  part->f.f[1] = 0;
  part->f.f[2] = 0;

#ifdef ROTATION
  {
    double scale;
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;

    /* and rescale quaternion, so it is exactly of unit length */	
    scale = sqrt( SQR(part->r.quat[0]) + SQR(part->r.quat[1]) +
		  SQR(part->r.quat[2]) + SQR(part->r.quat[3]));
    part->r.quat[0]/= scale;
    part->r.quat[1]/= scale;
    part->r.quat[2]/= scale;
    part->r.quat[3]/= scale;
  }
#endif
}

void init_forces()
{
  Cell *cell;
  Particle *p;
  int np, c, i;

  /* initialize forces with thermostat forces
     set torque to zero for all and rescale quaternions
  */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->pList.part;
    np = cell->pList.n;
    for (i = 0; i < np; i++)
      init_local_particle_force(&p[i]);
  }

  /* initialize ghost forces with zero
     set torque to zero for all and rescale quaternions
  */
  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    p  = cell->pList.part;
    np = cell->pList.n;
    for (i = 0; i < np; i++)
      init_ghost_force(&p[i]);
  }

#ifdef CONSTRAINTS
  init_constraint_forces();
#endif
}

#endif
