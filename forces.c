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
#include "pairwise.h"
#include "forces.h"


/************************************************************/

void force_init()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_init:\n",this_node));
  FORCE_TRACE(fprintf(stderr,"%d: found %d particles types\n",
		      this_node,n_particle_types));
}

/** nonbonded and bonded force calculation using the verlet list */
void calculate_verlet_ia()
{
  Cell *cell;
  Particle *p, **pairs;
  Particle *p1, *p2;
  int k, i,j, m, n, o, np;
  double d[3], dist2, dist;
  IA_parameters *ia_params;

  /* force calculation loop. */
  INNER_CELLS_LOOP(m, n, o) {
    cell = CELL_PTR(m, n, o);
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

	ia_params = get_ia_param(p1->r.type,p2->r.type);
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

void init_forces()
{
  Particle *p, *part;
  int np, m, n, o, i;
  int is_ghost;

  /* initialize forces with thermostat forces and
     ghost forces with zero
     set torque to zero for all and rescale quaternions
  */
  CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    is_ghost = IS_GHOST_CELL(m, n, o);

    for (i = 0; i < np; i++) {
      part = &p[i];
      if (is_ghost) {
	/* ghost particle selection */
	part->f[0] = 0;
	part->f[1] = 0;
	part->f[2] = 0;
      }
      else {
	/* real particle selection */
	friction_thermo(part);
#ifdef EXTERNAL_FORCES   
	if(part->ext_flag == PARTICLE_EXT_FORCE) {
	  part->f[0] += part->ext_force[0];
	  part->f[1] += part->ext_force[1];
	  part->f[2] += part->ext_force[2];
	}
#endif
      }

#ifdef ROTATION
      /* rotation for both ghost and real */
      {
	double scale;
	/* set torque to zero */
	part->torque[0] = 0;
	part->torque[1] = 0;
	part->torque[2] = 0;

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
  }

#ifdef CONSTRAINTS
  init_constraint_forces();
#endif
}
