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
#ifndef _LJGEN_H
#define _LJGEN_H

/** \file ljgen.hpp Routines to calculate the generalized lennard jones
 *  energy and/or force for a particle pair. "Generalized" here means
 *  that the LJ energy is of the form
 *
 *  eps * [ b1 * (sigma/(r-r_offset))^a1 - b2 * (sigma/(r-r_offset))^a2 + shift]
 *
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"

#ifdef LENNARD_JONES_GENERIC

int ljgen_set_params(int part_type_a, int part_type_b,
		     double eps, double sig, double cut,
		     double shift, double offset,
		     int a1, int a2, double b1, double b2,
		     double cap_radius
#ifdef LJGEN_SOFTCORE
         , double lambda, double softrad
#endif
         );

/** Calculate lennard Jones force between particle p1 and p2 */
inline void add_ljgen_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				   double d[3], double dist, double force[3])
{
  if (CUTOFF_CHECK(dist < ia_params->LJGEN_cut+ia_params->LJGEN_offset)) {
    int j;
    double r_off, frac, fac=0.0;
    r_off = dist - ia_params->LJGEN_offset;
#ifdef LJGEN_SOFTCORE
    r_off *= r_off;
    r_off += pow(ia_params->LJGEN_sig,2) * (1.0-ia_params->LJGEN_lambda)
             * ia_params->LJGEN_softrad;
    /* Taking a square root is not optimal, but we can't prevent the user from
       using an odd m, n coefficient. */
    r_off  = sqrt(r_off);
#endif
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJGEN_capradius) {
      frac = ia_params->LJGEN_sig/r_off;
      fac   = ia_params->LJGEN_eps 
#ifdef LJGEN_SOFTCORE
        * ia_params->LJGEN_lambda
#endif
	      * (ia_params->LJGEN_b1 * ia_params->LJGEN_a1 * pow(frac, ia_params->LJGEN_a1) 
	   - ia_params->LJGEN_b2 * ia_params->LJGEN_a2 * pow(frac, ia_params->LJGEN_a2))
	/ (r_off * dist);
      for(j=0;j<3;j++)
	force[j] += fac * d[j];
      
#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: LJ-Gen-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
    }
    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac = ia_params->LJGEN_sig/ia_params->LJGEN_capradius;
      fac   = ia_params->LJGEN_eps 
	* (ia_params->LJGEN_b1 * ia_params->LJGEN_a1 * pow(frac, ia_params->LJGEN_a1) 
	   - ia_params->LJGEN_b2 * ia_params->LJGEN_a2 * pow(frac, ia_params->LJGEN_a2))
	/ (ia_params->LJGEN_capradius * dist);
      for(j=0;j<3;j++)
	/* vector d is rescaled to length LJGEN_capradius */
	force[j] += fac * d[j];
    }
    /* this should not happen! */
    else {
      LJ_TRACE(fprintf(stderr, "%d: Lennard-Jones-Generic warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));
      
      frac = ia_params->LJGEN_sig/ia_params->LJGEN_capradius;
      fac   = ia_params->LJGEN_eps 
#ifdef LJGEN_SOFTCORE
        * ia_params->LJGEN_lambda
#endif
        * (ia_params->LJGEN_b1 * ia_params->LJGEN_a1 * pow(frac, ia_params->LJGEN_a1) 
	   - ia_params->LJGEN_b2 * ia_params->LJGEN_a2 * pow(frac, ia_params->LJGEN_a2))
	/ ia_params->LJGEN_capradius;
      
      force[0] += fac * ia_params->LJGEN_capradius;
    }
    
    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: LJGEN   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: LJGEN   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
    
    LJ_TRACE(fprintf(stderr,"%d: LJGEN: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
} 

/** calculate Lennard jones energy between particle p1 and p2. */
inline double ljgen_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac;

  if (CUTOFF_CHECK(dist < ia_params->LJGEN_cut+ia_params->LJGEN_offset)) {
    r_off = dist - ia_params->LJGEN_offset;
#ifdef LJGEN_SOFTCORE
    r_off *= r_off;
    r_off += pow(ia_params->LJGEN_sig,2) * (1.0-ia_params->LJGEN_lambda)
             * ia_params->LJGEN_softrad;
    /* Taking a square root is not optimal, but we can't prevent the user from
       using an odd m, n coefficient. */
    r_off  = sqrt(r_off);
#endif
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJGEN_capradius) {
      frac = ia_params->LJGEN_sig/r_off;
      return ia_params->LJGEN_eps
#ifdef LJGEN_SOFTCORE
        * ia_params->LJGEN_lambda
#endif
        *(ia_params->LJGEN_b1 * pow(frac, ia_params->LJGEN_a1) 
	      - ia_params->LJGEN_b2 * pow(frac, ia_params->LJGEN_a2)
				    + ia_params->LJGEN_shift);
    }
    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac = ia_params->LJGEN_sig/ia_params->LJGEN_capradius;
      return ia_params->LJGEN_eps
#ifdef LJGEN_SOFTCORE
        * ia_params->LJGEN_lambda
#endif
        *( ia_params->LJGEN_b1 * pow(frac, ia_params->LJGEN_a1) 
	      - ia_params->LJGEN_b2 * pow(frac, ia_params->LJGEN_a2)
				    + ia_params->LJGEN_shift);
    }
    /* this should not happen! */
    else {
      frac = ia_params->LJGEN_sig/ia_params->LJGEN_capradius;
      return ia_params->LJGEN_eps*(
	      ia_params->LJGEN_b1 * pow(frac, ia_params->LJGEN_a1) 
	       - ia_params->LJGEN_b2 * pow(frac, ia_params->LJGEN_a2)
				    + ia_params->LJGEN_shift);
    }
  }
  return 0.0;
}

/** calculate lj_capradius from force_cap */
void calc_ljgen_cap_radii();

#endif

/* LJGEN_H */
#endif 
