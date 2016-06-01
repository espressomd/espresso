
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

#ifdef DRUDE

/*
double rnd_kick_fac_d;
double rnd_kick_fac_c;
double fac_langevin_c_vc;
double fac_langevin_c_vd_and_d_vc;
double fac_langevin_d_vd;
*/

// set the parameters for the Drude correction
int drude_set_params(int bond_type, double temp_core, double gamma_core, double temp_drude, double gamma_drude, double k, double mass_drude, double r_cut);

// pre- or recalculate parameters for langevin cross terms
//void drude_recalc_params();


/** Computes the necessary corrections for the properly simulation of Drude Particles
    and adds this force to the particle forces (see \ref tclcommand_inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  Parameters of interaction
    @param dx        change in position
    @param forc, force1 and force2     force on particles
    @return true if bond is broken
*/

inline int calc_drude_forces(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force1[3], double force2[3])
{
  int dummy,i;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((iaparams->p.drude.r_cut > 0.0) && (dist > iaparams->p.drude.r_cut))
    return 1;

  double force_harmonic[3] = {0., 0., 0.};
  double force_subt_elec[3] = {0., 0., 0.};
  double force_lv_com[3] = {0., 0., 0.};
  double force_lv_dist[3] = {0., 0., 0.};

  double chgfac = p1->p.q*p2->p.q;

  double fac_harmonic = -iaparams->p.drude.k;
  double gamma_c = iaparams->p.drude.gamma_core;
  double gamma_d = iaparams->p.drude.gamma_drude;
  double temp_c = iaparams->p.drude.temp_core;
  double temp_d = iaparams->p.drude.temp_drude;

  double mass_c = p1->p.mass;
  double mass_d = p2->p.mass;
  double mass_tot = mass_d + mass_c;
  double mass_tot_inv = 1.0 / mass_tot;
  double mass_red = mass_d * mass_c / mass_tot;
  double mass_red_inv = 1.0 / mass_red;

  //double rnd_c[3] = { (d_random()-0.5), (d_random()-0.5), (d_random()-0.5) };
  //double rnd_d[3] = { (d_random()-0.5), (d_random()-0.5), (d_random()-0.5) };


  for (i=0;i<3;i++)  {
    double com_vel = mass_tot_inv * (mass_c * p1->m.v[i] + mass_d * p2->m.v[i]);
    //force_lv_com[i] =  -gamma_c / time_step * com_vel + sqrt(24.0 * gamma_c / time_step * temp_c) * (d_random()-0.5);
    force_lv_com[i] =  -gamma_c / time_step * com_vel + sqrt(2.0 * gamma_c / time_step * temp_c) * gaussian_random();
    double dist_vel = p2->m.v[i] - p1->m.v[i];
    //force_lv_dist[i] =  -gamma_d / time_step * dist_vel + sqrt(24.0 * gamma_d / time_step * temp_d) * (d_random()-0.5);
    force_lv_dist[i] =  -gamma_d / time_step * dist_vel + sqrt(2.0 * gamma_d / time_step * temp_d) * gaussian_random();
  }

  /* Apply forces: 
     -Harmonic bond
     -Langevin thermostat on distance core-drude and com in lab coords result in cross terms for velocities and rnd kicks */


  if (dist<ROUND_ERROR_PREC) {  /* dx[] == 0: the force is undefined. Let's use a random direction and no spring */
    for(i=0;i<3;i++) {
    	dx[i] = d_random()-0.5;
    }
  	dist2 = sqrlen(dx);
  	dist = sqrt(dist2);
    fac_harmonic = 0;
    fprintf(stderr,"dist<ROUND_ERROR_PREC");
  }
  
  for (i=0;i<3;i++)  {
     force_harmonic[i] = fac_harmonic*dx[i];
    
     force1[i] = mass_c * mass_tot_inv * force_lv_com[i] - force_lv_dist[i] + force_harmonic[i]; //Core
     force2[i] = mass_d * mass_tot_inv * force_lv_com[i] + force_lv_dist[i] - force_harmonic[i]; //Drude
     
     //force_subt_elec[i] = -coulomb.prefactor * chgfac * dx[i] / dist / dist2;
     //force1[i] = mass_c * mass_tot_inv * force_lv_com[i] - force_lv_dist[i] + force_harmonic[i] + force_subt_elec[i]; //Core
     //force2[i] = mass_d * mass_tot_inv * force_lv_com[i] + force_lv_dist[i] - force_harmonic[i] - force_subt_elec[i]; //Drude
     
     //force1[i] = force_com[i] + force_harmonic[i] + force_subt_elec[i]; //Core
     //force2[i] = force_com[i] - force_harmonic[i] - force_subt_elec[i]; //Drude
  }
  
  //fprintf(stderr,"Core Tot: %g %g %g\n", force1[0],force1[1],force1[2]);
  //fprintf(stderr,"Drude Tot: %g %g %g\n", force2[0],force2[1],force2[2]);

 /* 
  fprintf(stderr,"\ndx: %g %g %g\n", dx[0],dx[1],dx[2]);
  fprintf(stderr,"dv: %g %g %g\n", p2->m.v[0] - p1->m.v[0],p2->m.v[1] - p1->m.v[1],p2->m.v[2] - p1->m.v[2]);
  fprintf(stderr,"Harmonic: %g %g %g\n", force_harmonic[0],force_harmonic[1],force_harmonic[2]);
  fprintf(stderr,"Subt_elec: %g %g %g\n", force_subt_elec[0],force_subt_elec[1],force_subt_elec[2]);
  fprintf(stderr,"Dist: %g %g %g\n", force_lv_dist[0],force_lv_dist[1],force_lv_dist[2]);
  fprintf(stderr,"Com: %g %g %g\n", force_lv_com[0],force_lv_com[1],force_lv_com[2]);
  fprintf(stderr,"Core Tot: %g %g %g\n", force1[0],force1[1],force1[2]);
  fprintf(stderr,"Drude Tot: %g %g %g\n", force2[0],force2[1],force2[2]);
 */

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: DRUDE f = (%.3e,%.3e,%.3e)\n",this_node,p1->f.f[0]+force1[0],p1->f.f[1]+force1[1],p1->f.f[2]+force1[2]));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: DRUDE f = (%.3e,%.3e,%.3e)\n",this_node,p2->f.f[0]+force2[0],p2->f.f[1]+force2[1],p2->f.f[2]+force2[2]));
  return 0;
  
}


inline int drude_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double chgfac = p1->p.q*p2->p.q;

  if ((iaparams->p.drude.r_cut > 0.0) && 
      (dist > iaparams->p.drude.r_cut)) 
    return 1;
   //Harmonic 
  *_energy = 0.5*iaparams->p.drude.k*dist2;
   //Harmonic - Coulomb 
  //*_energy = 0.5*iaparams->p.drude.k*dist2 - coulomb.prefactor * chgfac / dist;
  
  return 0;
}

#endif
#endif
