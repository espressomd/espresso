/** \file forces.c Force calculation.
 *
 *  For more information see \ref forces.h "forces.h".
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "forces.h"
#include "debug.h"
#include "thermostat.h"
#include "communication.h"
#include "ghosts.h" 
#include "verlet.h"
#include "utils.h"
#include "interaction_data.h"
#include "grid.h"
#include "p3m.h"
#include "cells.h"
/* include the force files */
#include "fene.h"
#include "angle.h"

double minimum_part_dist = -1;

/************************************************************/

void force_init()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_init:\n",this_node));
  FORCE_TRACE(fprintf(stderr,"%d: found %d interaction types\n",
		      this_node,n_interaction_types));
  FORCE_TRACE(fprintf(stderr,"%d: found %d particles types\n",
		      this_node,n_particle_types));
}


/************************************************************/

static Particle *checked_particle_ptr(int id)
{
  Particle *p = local_particles[id];
  if(!p) {
    fprintf(stderr,"%d: ERROR: Atom %d has bond to unknown particle "
	    "(probably on different node)\n",this_node, id); 
    errexit();
  }
  return p;
}

/************************************************************/

void force_calc()
{
  Cell *cell;
  Particle *p, **pairs;
  Particle *p1, *p2;
  int k, i,j, m, n, o, np;
  double d[3], dist2, dist;
  IA_parameters *ia_params;
  double frac2,frac6,r_off;
  double fac, adist;
  /* bonded interactions */
  int type_num;
  /* electrostatic */
  double  erfc_part_ri;

  FORCE_TRACE(fprintf(stderr,"%d: force_calc: for %d (P %d,G %d)\n",this_node,n_particles+n_ghosts,n_particles,n_ghosts));

  minimum_part_dist = box_l[0] + box_l[1] + box_l[2];

  /* initialize forces with thermostat forces and
     ghost forces with zero */

  CELLS_LOOP(m, n, o) {
    p  = CELL_PTR(m, n, o)->pList.part;
    np = CELL_PTR(m, n, o)->pList.n;
    /* ghost selection */
    if (m == 0 || m == ghost_cell_grid[0] - 1 ||
	n == 0 || n == ghost_cell_grid[1] - 1 ||
	o == 0 || o == ghost_cell_grid[2] - 1) {
      for (i = 0; i < np; i++) {
	p[i].f[0] = 0;
	p[i].f[1] = 0;
	p[i].f[2] = 0;
      }
    }
    else {
      for (i = 0; i < np; i++)
	friction_thermo(&p[i]);
    }
  }

  INNER_CELLS_LOOP(m, n, o) {
    cell = CELL_PTR(m, n, o);
    p  = cell->pList.part;
    np = cell->pList.n;

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
	case BONDED_IA_ANGLE:
	  add_angle_pair_force(p1,
			       checked_particle_ptr(p1->bl.e[i+1]),
			       checked_particle_ptr(p1->bl.e[i+2]), type_num);
	  i+=3; break;
	default :
	  fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",p1->r.identity);
	  i = p1->bl.n; 
	  break;
	}
      }
    }

    /* calculate non bonded interactions (loop verlet lists of neighbors) */
    for (k = 0; k < cell->n_neighbors; k++) {
      pairs = cell->nList[k].vList.pair;
      np    = cell->nList[k].vList.n;
      for (i = 0; i < np; i++) {
	p1 = pairs[2*i];
	p2 = pairs[2*i+1];
	ia_params = get_ia_param(p1->r.type,p2->r.type);
	for(j=0;j<3;j++)
	  d[j] = p1->r.p[j] - p2->r.p[j];
	dist2 = SQR(d[0]) + SQR(d[1]) + SQR(d[2]);
	dist = sqrt(dist2);

	/* lennnard jones */
	if(dist < ia_params->LJ_cut+ia_params->LJ_offset) {
	  r_off = dist - ia_params->LJ_offset;
	  if(r_off>0.0) {
#ifdef LJ_WARN_WHEN_CLOSE
	    if (r_off < 0.9*ia_params->LJ_sig) {
	      fprintf(stderr, "Lennard-Jones warning: particles getting close\n");
	    }
#endif
	    frac2 = SQR(ia_params->LJ_sig/r_off);
	    frac6 = frac2*frac2*frac2;
	    fac = 48.* ia_params->LJ_eps * frac6*(frac6 - 0.5)*frac2;
	    for(j=0;j<3;j++) {
	      p1->f[j] += fac * d[j];
	      p2->f[j] -= fac * d[j];
	    }
	  }
	}

	/* real space coulomb */
	if(dist < p3m.r_cut) {
	  adist = p3m.alpha * dist;
	  erfc_part_ri = AS_erfc_part(adist) / dist;
	  fac = p3m.bjerrum * p1->r.q * p2->r.q  * 
	    exp(-adist*adist) * (erfc_part_ri + 2.0*p3m.alpha/1.772453851) / dist2;
	  p1->f[0] += fac * d[0];
	  p1->f[1] += fac * d[1];
	  p1->f[2] += fac * d[2];
	  p2->f[0] -= fac * d[0];
	  p2->f[1] -= fac * d[1];
	  p2->f[2] -= fac * d[2];
	}

	/* ramp */
	if(dist < ia_params->ramp_cut) {
	  if (dist < 1e-4) {
	    p1->f[0] += ia_params->ramp_force;
	    p1->f[1] += 0;
	    p1->f[2] += 0;
	    p2->f[0] -= ia_params->ramp_force;
	    p2->f[1] -= 0;
	    p2->f[2] -= 0;
	  }
	  else {
	    fac = ia_params->ramp_force/dist;
	    p1->f[0] += fac * d[0];
	    p1->f[1] += fac * d[1];
	    p1->f[2] += fac * d[2];
	    p2->f[0] -= fac * d[0];
	    p2->f[1] -= fac * d[1];
	    p2->f[2] -= fac * d[2];
	  }
	}
	/* minimal particle distance calculation */
	if (dist < minimum_part_dist)
	  minimum_part_dist = dist;
      } 
    }
  }

  /* calculate k-space part of electrostatic interaction. */
  if(p3m.bjerrum != 0.0) P3M_calc_kspace_forces();
}

/************************************************************/

void force_exit()
{
  FORCE_TRACE(fprintf(stderr,"%d: force_exit:\n",this_node));
}


