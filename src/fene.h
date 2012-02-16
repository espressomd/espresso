/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
#ifndef FENE_H
#define FENE_H
/** \file fene.h
 *  Routines to calculate the FENE Energy or/and FENE force 
 *  for a particle pair.
 *  \ref forces.c
 */

#include "utils.h"
#include "interaction_data.h"
#include "particle_data.h"
#include "random.h"

/************************************************************/

/// set the parameters for the fene potential
int fene_set_params(int bond_type, double k, double drmax, double r0);

/** Computes the FENE pair force and adds this
    force to the particle forces (see \ref interaction_data.c). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref interaction_data.c).
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return true if the bond is broken
*/
MDINLINE int calc_fene_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  int i;
  double fac, dr, len2, len;
 
  len2 = sqrlen(dx);
  len = sqrt(len2);
  dr = len - iaparams->p.fene.r0;

  if(dr >= iaparams->p.fene.drmax)
    return 1;

  fac = -iaparams->p.fene.k * dr / ((1.0 - dr*dr*iaparams->p.fene.drmax2i));
  if (fabs(dr) > ROUND_ERROR_PREC) {
     if(len > ROUND_ERROR_PREC) {  /* Regular case */
	fac /= len ; 
     } else { /* dx[] == 0: the force is undefined. Let's use a random direction */
        for(i=0;i<3;i++) dx[i] = d_random()-0.5;
        fac /= sqrt(sqrlen(dx));
     }
  } else { 
    fac = 0.0;
  }
  
  FENE_TRACE(if(fac > 50) fprintf(stderr,"WARNING: FENE force factor between Pair (%d,%d) large: %f at distance %f\n", p1->p.identity,p2->p.identity,fac,sqrt(len2)) );
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(len2),fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: FENE f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(len2),fac));
  
  return 0;
}

MDINLINE int fene_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double energy, dr;

  /* compute bond stretching (r-r0) */
  dr = sqrt(sqrlen(dx))-iaparams->p.fene.r0;

  /* check bond stretching */
  if(dr >= iaparams->p.fene.drmax) {
    char *errtext = runtime_error(128 + 2*ES_INTEGER_SPACE);
    ERROR_SPRINTF(errtext,"{077 FENE bond broken between particles %d and %d} ", p1->p.identity, p2->p.identity); 
    return 1;
  }

  energy = -0.5*iaparams->p.fene.k*iaparams->p.fene.drmax2;
  energy *= log((1.0 - dr*dr*iaparams->p.fene.drmax2i));
  *_energy = energy;
  return 0;
}

#endif
