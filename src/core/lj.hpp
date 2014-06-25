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
#ifndef _LJ_H
#define _LJ_H

/** \file lj.hpp
 *  Routines to calculate the lennard jones energy and/or  force 
 *  for a particle pair.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "mol_cut.hpp"
#include "forcecap.hpp"

#ifdef LENNARD_JONES
int ljforcecap_set_params(double ljforcecap);

int lennard_jones_set_params(int part_type_a, int part_type_b,
				      double eps, double sig, double cut,
				      double shift, double offset,
				      double cap_radius, double min);

/** Calculate lennard Jones force between particle p1 and p2 */
inline void add_lj_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  int j;
  double r_off, frac2, frac6, fac=0.0 ;
#ifdef SHANCHEN
  double order;
  double SixtRootOfTwo=1.12246204830937;
#endif 
  if (CUTOFF_CHECK(dist < ia_params->LJ_cut+ia_params->LJ_offset) &&
	  CUTOFF_CHECK(dist > ia_params->LJ_min+ia_params->LJ_offset))
  {
    r_off = dist - ia_params->LJ_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJ_capradius) {
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (r_off * dist);
#ifdef SHANCHEN
      if(ia_params->affinity_on==1){
         if(LB_COMPONENTS==2){
            if (dist > SixtRootOfTwo * ia_params->LJ_sig && 
                p1->r.composition[0] * p1->r.composition[1] * p2->r.composition[0] * p2->r.composition[1] > 0 ) {
               order = 2 * ( (p1->r.composition[0] - p1->r.composition[1])/(p1->r.composition[0] +  p1->r.composition[1])  + 
                       (p2->r.composition[0] - p2->r.composition[1])/(p2->r.composition[0] +  p2->r.composition[1]));
/* TODO: order should be rescaled properly by the lattice spacing */
               fac *= ( (1-ia_params->affinity[1]) * 0.5 * (1+tanh(-order)) + (1-ia_params->affinity[0]) *  0.5 * (1+tanh(order)) )  ;
            }
         }
      }
#endif
      for(j=0;j<3;j++)
	force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: LJ-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
    }
    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / (ia_params->LJ_capradius * dist);
      for(j=0;j<3;j++)
	/* vector d is rescaled to length LJ_capradius */
	force[j] += fac * d[j];
    }
    /* this should not happen! */
    else {
      LJ_TRACE(fprintf(stderr, "%d: Lennard-Jones warning: Particles id1=%d id2=%d exactly on top of each other\n",this_node,p1->p.identity,p2->p.identity));

      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJ_eps * frac6*(frac6 - 0.5) / ia_params->LJ_capradius;

      force[0] += fac * ia_params->LJ_capradius;
    }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: LJ   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    LJ_TRACE(fprintf(stderr,"%d: LJ: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
}

/** calculate Lennard jones energy between particle p1 and p2. */
inline double lj_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac2, frac6;
  if (CUTOFF_CHECK(dist < ia_params->LJ_cut+ia_params->LJ_offset) &&
	  CUTOFF_CHECK(dist > ia_params->LJ_min+ia_params->LJ_offset))
  {
    r_off = dist - ia_params->LJ_offset;
    /* normal case: resulting force/energy smaller than capping. */
    if(r_off > ia_params->LJ_capradius) {
      frac2 = SQR(ia_params->LJ_sig/r_off);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
    /* capped part of lj potential. */
    else if(dist > 0.0) {
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
    /* this should not happen! */
    else {
      frac2 = SQR(ia_params->LJ_sig/ia_params->LJ_capradius);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJ_eps*(SQR(frac6)-frac6+ia_params->LJ_shift);
    }
  }
  return 0.0;
}

/** calculate lj_capradius from force_cap */
void calc_lj_cap_radii();

#endif /* ifdef LENNARD_JONES */
#endif

