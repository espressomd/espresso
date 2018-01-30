
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

#ifndef DRUDE_H
#define DRUDE_H
/** \file drude_correction.h
 *  Routines to subtract the electrostatic Energy and/or the electrostatic force 
 *  for a particle pair and correct the equations of motion for a Drude pair.
 *  \ref forces.c
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "integrate.hpp"
#include "random.hpp"
#include "p3m.hpp"

#ifdef DRUDE

// Set the parameters for the Drude correction
int drude_set_params(int bond_type, double temp_com, double gamma_com, double temp_drude, double gamma_drude, double k, double r_cut);

/** Computes the necessary corrections for the properly simulation of Drude Particles
    and adds this force to the particle forces. 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  Parameters of interaction
    @param dx        change in position
    @param force1 and force2     force on particles
    @return true if bond is broken
*/

inline int calc_drude_forces(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force1[3], double force2[3])
{
  //Bond broke?
  double dist = Utils::veclen(dx);
  if (iaparams->p.drude.r_cut > 0.0 && dist > iaparams->p.drude.r_cut) {
    return 1;
  }

  double force_harmonic, force_lv_com, force_lv_dist, com_vel, dist_vel;
  double chgfac = p1->p.q*p2->p.q;
  double fac_harmonic = -iaparams->p.drude.k;
  double mass_tot_inv = 1.0 / (p1->p.mass + p2->p.mass);

  double force_p3m_sr[3] = {0,0,0}; 

  p3m_add_pair_force(chgfac, dx, dist*dist, dist, force_p3m_sr);

  for (int i=0; i<3; i++) {

      //Langevin thermostat for center of mass
      com_vel = mass_tot_inv * (p1->p.mass * p1->m.v[i] + p2->p.mass * p2->m.v[i]);
      //force_lv_com =  -iaparams->p.drude.gamma_com / time_step * com_vel + sqrt(2.0 * iaparams->p.drude.gamma_com / time_step * iaparams->p.drude.temp_com) * gaussian_random_cut();
      force_lv_com =  -iaparams->p.drude.gamma_com / time_step * com_vel + sqrt(24.0 * iaparams->p.drude.gamma_com / time_step * iaparams->p.drude.temp_com) * (d_random()-0.5);
  
      //Langevin thermostat for distance core->drude
      dist_vel = p2->m.v[i] - p1->m.v[i];
      //force_lv_dist =  -iaparams->p.drude.gamma_drude / time_step * dist_vel + sqrt(2.0 * iaparams->p.drude.gamma_drude / time_step * iaparams->p.drude.temp_drude) * gaussian_random_cut();
      force_lv_dist =  -iaparams->p.drude.gamma_drude / time_step * dist_vel + sqrt(24.0 * iaparams->p.drude.gamma_drude / time_step * iaparams->p.drude.temp_drude) * (d_random()-0.5);

      //Spring
      force_harmonic = fac_harmonic * dx[i];

      //Add forces
      force1[i] = p1->p.mass * mass_tot_inv * force_lv_com - force_lv_dist + force_harmonic - force_p3m_sr[i]; //Core
      force2[i] = p2->p.mass * mass_tot_inv * force_lv_com + force_lv_dist - force_harmonic + force_p3m_sr[i]; //Drude

  }

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: DRUDE f = (%.3e,%.3e,%.3e)\n",this_node,p1->f.f[0]+force1[0],p1->f.f[1]+force1[1],p1->f.f[2]+force1[2]));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: DRUDE f = (%.3e,%.3e,%.3e)\n",this_node,p2->f.f[0]+force2[0],p2->f.f[1]+force2[1],p2->f.f[2]+force2[2]));
  return 0;
  
}


inline int drude_energy(const Particle *p1, const Particle *p2, const Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double chgfac = p1->p.q*p2->p.q;

  if (iaparams->p.drude.r_cut > 0.0 && dist > iaparams->p.drude.r_cut) 
    return 1;
   //Harmonic, subtract drude-core shortrange energy
  *_energy = 0.5*iaparams->p.drude.k*dist2 + p3m_pair_energy(-chgfac, dist);
  
  return 0;
}

#endif
#endif
