/** \file forces.c Force calculation.
 *
 *  For more information see \ref forces.h "forces.h".
*/
#include <stdio.h>
#include <stdlib.h>
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


/************************************************************/

void force_init()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_init:\n",this_node));
  FORCE_TRACE(fprintf(stderr,"%d: found %d particles types\n",
		      this_node,n_particle_types));
}


void force_calc()
{
  Cell *cell;
  Particle *p, **pairs;
  Particle *p1, *p2;
  int k, i,j, m, n, o, np;
  double d[3], dist2, dist;
  IA_parameters *ia_params;
  /* preparation */
  init_forces();
#ifdef ROTATION
  init_torques();
#endif


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

  calc_long_range_forces();
}

/************************************************************/

void calc_bonded_forces(Particle *p, int np)
{
  Particle *p1;
  int i, j, type_num;

  /* calculate bonded interactions (loop local particles) */
  for(j = 0; j < np; j++) {
    p1 = &p[j];
    i=0;
    while(i<p1->bl.n) {
      type_num = p1->bl.e[i];
      switch(bonded_ia_params[type_num].type) {
      case BONDED_IA_FENE:
	add_fene_pair_force(p1,
			    checked_particle_ptr(p1->bl.e[i+1]), type_num);
	i+=2; break;
      case BONDED_IA_HARMONIC:
	add_harmonic_pair_force(p1,
				checked_particle_ptr(p1->bl.e[i+1]), type_num);
	i+=2; break;
      case BONDED_IA_ANGLE:
	add_angle_force(p1,
			checked_particle_ptr(p1->bl.e[i+1]),
			checked_particle_ptr(p1->bl.e[i+2]), type_num);
	i+=3; break;
      default :
	fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",p1->r.identity);
	i = p1->bl.n; 
	break;
      }
    }
#ifdef CONSTRAINTS
    add_constraints_forces(p1);
#endif
  }
}

/************************************************************/

void calc_long_range_forces()
{
#ifdef ELECTROSTATICS  
  /* calculate k-space part of electrostatic interaction. */
  if(coulomb.method == COULOMB_P3M) 
    P3M_calc_kspace_forces(1,0);
  else if (coulomb.method == COULOMB_MMM1D)
    MMM1D_calc_forces();
#endif
}

/************************************************************/

void init_forces()
{
  Particle *p;
  int np, m, n, o, i;
  /* initialize forces with thermostat forces and
     ghost forces with zero */
  CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    /* ghost particle selection */
    if (IS_GHOST_CELL(m, n, o)) {
      for (i = 0; i < np; i++) {
	p[i].f[0] = 0;
	p[i].f[1] = 0;
	p[i].f[2] = 0;
      }
    }
    /* real particle selection */
    else {
      for (i = 0; i < np; i++) {
       	friction_thermo(&p[i]);
#ifdef EXTERNAL_FORCES   
	if(p[i].ext_flag == PARTICLE_EXT_FORCE) {
	  p[i].f[0] += p[i].ext_force[0];
	  p[i].f[1] += p[i].ext_force[1];
	  p[i].f[2] += p[i].ext_force[2];
	}
#endif
      }
    }
  }
#ifdef CONSTRAINTS
  init_constraint_forces();
#endif
}


