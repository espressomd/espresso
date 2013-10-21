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
#ifndef _MORSE_H
#define _MORSE_H

/** \file morse.hpp
 *  Routines to calculate the lennard jones energy and/or  force 
 *  for a particle pair.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"
#include "forcecap.hpp"

#ifdef MORSE

int morse_set_params(int part_type_a, int part_type_b,
		     double eps, double alpha,
		     double rmin, double cut, double cap_radius);

/** calculate morse_capradius from force_cap */
void calc_morse_cap_radii();

/** Calculate Morse force between particle p1 and p2 */
inline void add_morse_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				   double d[3], double dist, double force[3])
{
  double add1, add2, fac=0.0;
  int j;
  if(CUTOFF_CHECK(dist < ia_params->MORSE_cut)) { 
    /* normal case: resulting force/energy smaller than capping. */

    if(dist > ia_params->MORSE_capradius ) {
      add1 = exp(-2.0 * ia_params->MORSE_alpha * (dist - ia_params->MORSE_rmin));
      add2 = exp( -ia_params->MORSE_alpha * (dist - ia_params->MORSE_rmin));
      fac   = -ia_params->MORSE_eps * 2.0 * ia_params->MORSE_alpha * (add2 - add1) / dist;

      for(j=0;j<3;j++)
	force[j] += fac * d[j];
#ifdef MORSE_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: Morse-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif

    }
    /* capped part of morse potential. */

    else if(dist > 0.0) {

        add1 = exp(-2.0 * ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        add2 = exp( -ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        fac   = -ia_params->MORSE_eps * 2.0 * ia_params->MORSE_alpha * (add2 - add1) / ia_params->MORSE_capradius;

	/* vector d is rescaled to length MORSE_capradius */
        for(j=0;j<3;j++)
	  force[j] += fac * d[j];
    }
    /* this should not happen! */

    else {
      MORSE_TRACE(fprintf(stderr, "%d: Morse warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));

      add1 = exp(-2.0 * ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
      add2 = exp( -ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
      fac   = -ia_params->MORSE_eps * 2.0 * ia_params->MORSE_alpha * (add2 - add1) / ia_params->MORSE_capradius;
      
      force[0] += fac * ia_params->MORSE_capradius;
    }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: MORSE   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: MORSE   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    MORSE_TRACE(fprintf(stderr,"%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));

  }
}

/** calculate Morse energy between particle p1 and p2. */
inline double morse_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{

  double add1, add2, fac;

  if(CUTOFF_CHECK(dist < ia_params->MORSE_cut)) {

    /* normal case: resulting force/energy smaller than capping. */
    if(dist > ia_params->MORSE_capradius) { 
        add1 = exp(-2.0 * ia_params->MORSE_alpha * (dist - ia_params->MORSE_rmin));
        add2 = 2.0 * exp( -ia_params->MORSE_alpha * (dist - ia_params->MORSE_rmin));
        fac   = ia_params->MORSE_eps * (add1 - add2) - ia_params->MORSE_rest;
      return fac;
    }

    /* capped part of morse potential. */
    else if(dist > 0.0) {
        add1 = exp(-2.0 * ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        add2 = 2.0 * exp( -ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        fac   = ia_params->MORSE_eps * (add1 - add2) - ia_params->MORSE_rest;
      return fac;
    }
    /* this should not happen! */

    else {
        add1 = exp(-2.0 * ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        add2 = 2.0 * exp( -ia_params->MORSE_alpha * (ia_params->MORSE_capradius - ia_params->MORSE_rmin));
        fac   = ia_params->MORSE_eps * (add1 - add2) - ia_params->MORSE_rest;
        return fac;
    }

  }

  return 0.0;
}

#endif /* ifdef MORSE */
#endif /* ifdef _MORSE_H */
