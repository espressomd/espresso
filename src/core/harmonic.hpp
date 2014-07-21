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
#ifndef _HARMONIC_HPP
#define _HARMONIC_HPP
/** \file harmonic.hpp
 *  Routines to calculate the HARMONIC Energy or/and HARMONIC force 
 *  for a particle pair.
 *  \ref forces.cpp
*/

/************************************************************/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "random.hpp"

/// set the parameters for the harmonic potential
int harmonic_set_params(int bond_type, double k, double r,double r_cut);

/** Computes the HARMONIC pair force and adds this
    force to the particle forces (see \ref interaction_data.cpp). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond type number of the angle interaction (see \ref interaction_data.cpp).
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return 0.
*/
inline int calc_harmonic_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  if ((iaparams->p.harmonic.r_cut > 0.0) &&
      (dist > iaparams->p.harmonic.r_cut)) 
    return 1;

  dr = dist - iaparams->p.harmonic.r;
  fac = -iaparams->p.harmonic.k * dr;
  if (fabs(dr) > ROUND_ERROR_PREC) {
     if(dist>ROUND_ERROR_PREC) {  /* Regular case */
        fac /= dist;
     } else { /* dx[] == 0: the force is undefined. Let's use a random direction */
        for(i=0;i<3;i++) dx[i] = d_random()-0.5;
        fac /= sqrt(sqrlen(dx));
     }
  } else { 
     fac=0;
  }
  
  for(i=0;i<3;i++)
    force[i] = fac*dx[i];
  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist2,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: HARMONIC f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist2,fac));

  return 0;
}

inline int harmonic_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((iaparams->p.harmonic.r_cut > 0.0) && 
      (dist > iaparams->p.harmonic.r_cut)) 
    return 1;

  *_energy = 0.5*iaparams->p.harmonic.k*SQR(dist - iaparams->p.harmonic.r);
  return 0;
}

#endif
