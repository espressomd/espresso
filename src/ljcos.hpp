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
#ifndef _LJCOS_H
#define _LJCOS_H
/** \file ljcos.hpp
 *  Routines to calculate the lennard jones+cosine energy and/or force 
 *  for a particle pair.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"

#ifdef LJCOS

int lj_cos_set_params(int part_type_a, int part_type_b,
		      double eps, double sig, double cut,
		      double offset);

inline void add_ljcos_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				   double d[3], double dist, double force[3])
{
  int j;
  double r_off, frac2, frac6, fac=0.0;

  if(CUTOFF_CHECK(dist < ia_params->LJCOS_cut+ia_params->LJCOS_offset)) {
    r_off = dist - ia_params->LJCOS_offset;
    /* cos part of ljcos potential. */
    if(dist > ia_params->LJCOS_rmin+ia_params->LJCOS_offset) {
      fac   = (r_off/dist) * ia_params->LJCOS_alfa * ia_params->LJCOS_eps * (sin(ia_params->LJCOS_alfa * SQR(r_off) + ia_params->LJCOS_beta));
      for(j=0;j<3;j++)
	force[j] += fac * d[j];
    }
    /* lennard-jones part of the potential. */
    else if(dist > 0) {
      frac2 = SQR(ia_params->LJCOS_sig/r_off);
      frac6 = frac2*frac2*frac2;
      fac   = 48.0 * ia_params->LJCOS_eps * frac6*(frac6 - 0.5) / (r_off * dist);

      for(j=0;j<3;j++)
	force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: LJCOS-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
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


inline double ljcos_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off, frac2, frac6;

  if(CUTOFF_CHECK(dist < ia_params->LJCOS_cut+ia_params->LJCOS_offset)) {
    r_off = dist-ia_params->LJCOS_offset;
    /* lennard-jones part of the potential. */
    if (dist < (ia_params->LJCOS_rmin+ia_params->LJCOS_offset)) {
      //printf("this is nomal ,  %.3e \n",r_off);
      frac2 = SQR(ia_params->LJCOS_sig/r_off);
      frac6 = frac2*frac2*frac2;
      return 4.0*ia_params->LJCOS_eps*(SQR(frac6)-frac6);
    }
    /* cosine part of the potential. */
    else if (dist < (ia_params->LJCOS_cut+ia_params->LJCOS_offset)) {
      return .5*ia_params->LJCOS_eps*(cos(ia_params->LJCOS_alfa*SQR(r_off)+ia_params->LJCOS_beta)-1.);
    }
    /* this should not happen! */
    else {
      fprintf(stderr,"this is the distance, which is negative %.3e\n",r_off);
    }
  }
  return 0.0;
}

#endif
#endif
