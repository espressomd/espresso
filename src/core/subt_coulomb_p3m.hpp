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
#ifndef SUBT_COULOMB_P3M_H
#define SUBT_COULOMB_P3M_H
/** \file subt_coublom_p3m.hpp
 *  Routines to subtract the Coulomb Energy and/or the Coulomb force 
 *  for a particle pair. Uses the exact short-ranged part of p3m and 
 *  precalculated Ewald sum for the long range part. As it is a bond,
 *  it can only work for up to max_r_cut. 
 *  \ref forces.cpp
 */

#include "config.hpp"

#ifdef P3M

#include "debug.hpp"
#include "utils.hpp"
#include "interaction_data.hpp"
#include "p3m.hpp"
#include "grid.hpp"

/// set the parameters for the subtract_coulomb_p3m potential
int subt_coulomb_p3m_set_params(int bond_type, int res, double r_max);

/** Computes the negative of the LENNARD-JONES pair forces 
  and adds this force to the particle forces (see \ref tclcommand_inter). 
  @param p1        Pointer to first particle.
  @param p2        Pointer to second/middle particle.
  @param iaparams  Parameters of interaction
  @param dx        change in position
  @param force     force on particles
  @return true if bond is broken
 */
inline int calc_subt_coulomb_p3m_pair_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3])
{
    IA_parameters *ia_params;
    double dist2 = sqrlen(dx);
    double dist = sqrt(dist2);
    double r_max = iaparams->p.subt_coulomb_p3m.r_max;
    double res = iaparams->p.subt_coulomb_p3m.res;

    double bin_size = 1.0*r_max/res;
    if (dist >= r_max - bin_size)
        return 1;
    else if (dist > 0) {

        double force_subtr[3] = {0,0,0};
        double q1q2 = p1->p.q*p2->p.q;

        //Add short range to force
        if (q1q2)
            p3m_add_pair_force(q1q2, dx, dist2, dist, force_subtr);

        //Add long range from precalculated array
        int dist_bin = dist/r_max * res; //The lower bin index of the array
        double d_bin = r_max * dist_bin / res; //Position of the lower bin
        double sp = (dist - d_bin) / bin_size; //Splitting factor to next bin for linear interpol.
        //Long range calulation: e(p1,p2) * q1*q2/V * F_Ewald_precomputed
        double lr = q1q2 / (box_l[0]*box_l[1]*box_l[2]) / dist * ((1.0-sp) * iaparams->p.subt_coulomb_p3m.long_range_forces[dist_bin] + sp * iaparams->p.subt_coulomb_p3m.long_range_forces[dist_bin+1]);
        for(int i=0; i<3; i++)
        {
            force_subtr[i] += dx[i] * lr;
            //Invert
            force[i] -= force_subtr[i];
        }

        ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_COULOMB_P3M f = (%.3e,%.3e,%.3e) with part id=%d at dist %f\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,sqrt(dist2)));
        ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: SUBT_COULOMB_P3M f = (%.3e,%.3e,%.3e) with part id=%d at dist %f\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,sqrt(dist2)));
    }
    return 0;
}

inline int subt_coulomb_p3m_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
    double energy=0.0;
    IA_parameters *ia_params;
    //double dist2 = sqrlen(dx);
    //double dist = sqrt(dist2);

    //ia_params = get_ia_param(p1->p.type,p2->p.type);

    *_energy = energy;

    return 0;
}

#endif

#endif
