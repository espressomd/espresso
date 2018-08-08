/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _HARMONIC_DUMBBELL_HPP
#define _HARMONIC_DUMBBELL_HPP
/** \file harmonic_dumbbell.hpp
 *  Routines to calculate the HARMONIC Energy or/and HARMONIC force 
 *  for a particle pair.
 *  \ref forces.cpp
*/

/************************************************************/

#include "debug.hpp"
#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "random.hpp"

#ifdef ROTATION

/// set the parameters for the harmonic potential
int harmonic_dumbbell_set_params(int bond_type, double k1, double k2, double r, double r_cut);

/** Computes the HARMONIC pair force and adds this
    force to the particle forces (see \ref interaction_data.cpp). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref interaction_data.cpp).
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return 0.
*/
inline int calc_harmonic_dumbbell_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  if ((iaparams->p.harmonic_dumbbell.r_cut > 0.0) &&
      (dist > iaparams->p.harmonic_dumbbell.r_cut)) 
    return 1;

  dr = dist - iaparams->p.harmonic_dumbbell.r;
  fac = -iaparams->p.harmonic_dumbbell.k1 * dr;
  if (fabs(dr) > ROUND_ERROR_PREC) {
     if (dist > ROUND_ERROR_PREC)  /* Regular case */
        fac /= dist;
     else { /* dx[] == 0: the force is undefined. Let's use a random direction */
        for(int i=0;i<3;i++)
	      dx[i] = d_random()-0.5;
        fac /= sqrt(sqrlen(dx));
     }
  } else { 
     fac = 0;
  }
  
  for (int i=0; i<3; i++)
    force[i] = fac*dx[i];

  double dhat[3];
  dhat[0] = dx[0]/dist;
  dhat[1] = dx[1]/dist;
  dhat[2] = dx[2]/dist;

  double da[3];
  da[0] = dhat[1]*p1->r.quatu[2] - dhat[2]*p1->r.quatu[1];
  da[1] = dhat[2]*p1->r.quatu[0] - dhat[0]*p1->r.quatu[2];
  da[2] = dhat[0]*p1->r.quatu[1] - dhat[1]*p1->r.quatu[0];

  p1->f.torque[0] += iaparams->p.harmonic_dumbbell.k2 * da[0];
  p1->f.torque[1] += iaparams->p.harmonic_dumbbell.k2 * da[1];
  p1->f.torque[2] += iaparams->p.harmonic_dumbbell.k2 * da[2];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));

  return 0;
}

inline int harmonic_dumbbell_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((iaparams->p.harmonic_dumbbell.r_cut > 0.0) && 
      (dist > iaparams->p.harmonic_dumbbell.r_cut)) 
    return 1;

  double dhat[3];
  dhat[0] = dx[0]/dist;
  dhat[1] = dx[1]/dist;
  dhat[2] = dx[2]/dist;

  double da[3];
  da[0] = dhat[1]*p1->r.quatu[2] - dhat[2]*p1->r.quatu[1];
  da[1] = dhat[2]*p1->r.quatu[0] - dhat[0]*p1->r.quatu[2];
  da[2] = dhat[0]*p1->r.quatu[1] - dhat[1]*p1->r.quatu[0];

  double torque[3];
  torque[0] = iaparams->p.harmonic_dumbbell.k2 * da[0];
  torque[1] = iaparams->p.harmonic_dumbbell.k2 * da[1];
  torque[2] = iaparams->p.harmonic_dumbbell.k2 * da[2];

  double diff[3];
  diff[0] = dhat[0] - p1->r.quatu[0];
  diff[1] = dhat[1] - p1->r.quatu[1];
  diff[2] = dhat[2] - p1->r.quatu[2];

  *_energy = 0.5*iaparams->p.harmonic_dumbbell.k1*Utils::sqr(dist - iaparams->p.harmonic.r)
           + 0.5*iaparams->p.harmonic_dumbbell.k2*(torque[0]*diff[0] + torque[1]*diff[1] + torque[2]*diff[2]);
  return 0;
}

#endif // ROTATION

#endif
