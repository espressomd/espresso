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
#ifndef DPD_H
#define DPD_H
/** \file dpd.hpp
 *  Routines to use dpd as thermostat or pair force
 *  T. Soddemann, B. Duenweg and K. Kremer, Phys. Rev. E 68, 046702 (2003)
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "thermostat.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "virtual_sites.hpp"

/** Flag to decide wether to allow for fixed particles with DPD */
extern int dpd_ignore_fixed_particles;

/** DPD Friction coefficient gamma. */
extern double dpd_gamma;
/** DPD thermostat cutoff */
extern double dpd_r_cut;
/** DPD thermostat weight function */
extern int dpd_wf;

/** DPD transversal Friction coefficient gamma. */
extern double dpd_tgamma;
/** trans DPD thermostat cutoff */
extern double dpd_tr_cut;
/** trans DPD thermostat weight function */
extern int dpd_twf;

#ifdef DPD
extern double dpd_r_cut_inv;
extern double dpd_pref1;
extern double dpd_pref2;

#ifdef TRANS_DPD 
extern double dpd_tr_cut_inv;
extern double dpd_pref3;
extern double dpd_pref4;
#endif

void dpd_switch_off(void);
void thermo_init_dpd();
void dpd_heat_up();
void dpd_cool_down();

/** Calculate Random Force and Friction Force acting between particle
    p1 and p2 and add them to their forces. */
inline void add_dpd_thermo_pair_force(Particle *p1, Particle *p2, double d[3], double dist, double dist2)
{
  extern double dpd_gamma,dpd_pref1, dpd_pref2,dpd_r_cut,dpd_r_cut_inv;
  extern int dpd_wf;
#ifdef TRANS_DPD
  extern double dpd_tgamma, dpd_pref3, dpd_pref4,dpd_tr_cut,dpd_tr_cut_inv;
  extern int dpd_twf;
#endif
  int j;
  // velocity difference between p1 and p2
  double vel12_dot_d12=0.0;
  // inverse distance
  double dist_inv;
  // weighting functions for friction and random force
  double omega,omega2;// omega = w_R/dist
  double friction, noise;
  //Projection martix
#ifdef TRANS_DPD
  int i;
  double P_times_dist_sqr[3][3]={{dist2,0,0},{0,dist2,0},{0,0,dist2}},noise_vec[3];
  double f_D[3],f_R[3];
#endif
  double tmp;
#ifdef DPD_MASS
  double massf;
#endif

#ifdef EXTERNAL_FORCES
  // if any of the two particles is fixed in some direction then
  // do not add any dissipative or stochastic dpd force part
  // because dissipation-fluctuation theorem is violated
  if (dpd_ignore_fixed_particles)
    if ( (p1->l.ext_flag | p2->l.ext_flag) & COORDS_FIX_MASK) return;
#endif

#ifdef VIRTUAL_SITES
    if (ifParticleIsVirtual(p1) || ifParticleIsVirtual(p2)) return;
#endif	  

#ifdef DPD_MASS_RED
  massf=2*PMASS(*p1)*PMASS(*p2)/(PMASS(*p1)+PMASS(*p2));
#endif

#ifdef DPD_MASS_LIN
  massf=0.5*(PMASS(*p1)+PMASS(*p2));
#endif


  dist_inv = 1.0/dist;

  if((dist < dpd_r_cut)&&(dpd_gamma > 0.0)) {
    if ( dpd_wf == 1 ) //w_R=1
    {
       omega    = dist_inv;
    }
    else //w_R=(1-r/r_c)
    {
    	omega    = dist_inv- dpd_r_cut_inv;
    }
#ifdef DPD_MASS
    omega*=sqrt(massf);
#endif
    omega2   = SQR(omega);
    //DPD part
    // friction force prefactor
    for(j=0; j<3; j++)  vel12_dot_d12 += (p1->m.v[j] - p2->m.v[j]) * d[j];
    friction = dpd_pref1 * omega2 * vel12_dot_d12;
    // random force prefactor
    noise    = dpd_pref2 * omega      * (d_random()-0.5);
    for(j=0; j<3; j++) {
       p1->f.f[j] += ( tmp = (noise - friction)*d[j] );
       p2->f.f[j] -= tmp;
    }
  }
#ifdef TRANS_DPD
    //DPD2 part
  if ((dist < dpd_tr_cut)&&(dpd_tgamma > 0.0)){
      if ( dpd_twf == 1 )
      {
        omega    = dist_inv;
      }
      else 
      {
        omega    = dist_inv- dpd_tr_cut_inv;
      }
#ifdef DPD_MASS
      omega*=sqrt(massf);
#endif
      omega2   = SQR(omega);
      for (i=0;i<3;i++){
        //noise vector
        noise_vec[i]=d_random()-0.5;
        // Projection Matrix
        for (j=0;j<3;j++){
          P_times_dist_sqr[i][j]-=d[i]*d[j];
        }
      }
      for (i=0;i<3;i++){
        //Damping force
        f_D[i]=0;
        //Random force
        f_R[i]=0;
        for (j=0;j<3;j++){
          f_D[i]+=P_times_dist_sqr[i][j]*(p1->m.v[j] - p2->m.v[j]);
          f_R[i]+=P_times_dist_sqr[i][j]*noise_vec[j];
        }
        //NOTE: velocity are scaled with time_step
        f_D[i]*=dpd_pref3*omega2;
        //NOTE: noise force scales with 1/sqrt(time_step
        f_R[i]*=dpd_pref4*omega*dist_inv;
      }
      for(j=0; j<3; j++) {
        tmp=f_R[j]-f_D[j];
        p1->f.f[j] += tmp;
        p2->f.f[j] -= tmp;
      }
  }
#endif
}
#endif

#ifdef INTER_DPD
void inter_dpd_heat_up();
void inter_dpd_cool_down();
void inter_dpd_switch_off();
int inter_dpd_set_params(int part_type_a, int part_type_b,
			 double gamma, double r_c, int wf,
			 double tgamma, double tr_c,
			 int twf);
void inter_dpd_init();
void inter_dpd_update_params(double pref2_scale);

inline void add_inter_dpd_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double dist2)
{
  int j;
  // velocity difference between p1 and p2
  double vel12_dot_d12=0.0;
  // inverse distance
  double dist_inv;
  // weighting functions for friction and random force
  double omega,omega2;// omega = w_R/dist
  double friction, noise;
  //Projection martix
  int i;
  double P_times_dist_sqr[3][3]={{0,0,0},{0,0,0},{0,0,0}},noise_vec[3];
  double f_D[3],f_R[3];
  double tmp;
#ifdef DPD_MASS
  double massf;
#endif

#ifdef EXTERNAL_FORCES
  // if any of the two particles is fixed in some direction then
  // do not add any dissipative or stochastic dpd force part
  // because dissipation-fluctuation theorem is violated
  if (dpd_ignore_fixed_particles)
    if ( (p1->l.ext_flag | p2->l.ext_flag) & COORDS_FIX_MASK) return;
#endif

#ifdef DPD_MASS_RED
  massf=2*PMASS(*p1)*PMASS(*p2)/(PMASS(*p1)+PMASS(*p2));
#endif

#ifdef DPD_MASS_LIN
  massf=0.5*(PMASS(*p1)+PMASS(*p2));
#endif

  P_times_dist_sqr[0][0]=dist2;
  P_times_dist_sqr[1][1]=dist2;
  P_times_dist_sqr[2][2]=dist2;
  dist_inv = 1.0/dist;
  if((dist < ia_params->dpd_r_cut)&&(ia_params->dpd_gamma > 0.0)) {
    if ( dpd_wf == 1 )
    {
       omega    = dist_inv;
    }
    else 
    {
    	omega    = dist_inv - 1.0/ia_params->dpd_r_cut;
    }
#ifdef DPD_MASS
    omega*=sqrt(massf);
#endif
    omega2   = SQR(omega);
    //DPD part
       // friction force prefactor
    for(j=0; j<3; j++)  vel12_dot_d12 += (p1->m.v[j] - p2->m.v[j]) * d[j];
    friction = ia_params->dpd_pref1 * omega2 * vel12_dot_d12;
    // random force prefactor
    noise    = ia_params->dpd_pref2 * omega      * (d_random()-0.5);
    for(j=0; j<3; j++) {
       p1->f.f[j] += ( tmp = (noise - friction)*d[j] );
       p2->f.f[j] -= tmp;
    }
  }
  //DPD2 part
  if ((dist < ia_params->dpd_tr_cut)&&(ia_params->dpd_tgamma > 0.0)){
      if ( ia_params->dpd_twf == 1 )
      {
        omega    = dist_inv;
      }
      else 
      {
        omega    = dist_inv- 1.0/ia_params->dpd_tr_cut;
      }
#ifdef DPD_MASS
      omega*=sqrt(massf);
#endif
      omega2   = SQR(omega);
      for (i=0;i<3;i++){
        //noise vector
        noise_vec[i]=d_random()-0.5;
        // Projection Matrix
        for (j=0;j<3;j++){
          P_times_dist_sqr[i][j]-=d[i]*d[j];
        }
      }
      for (i=0;i<3;i++){
        //Damping force
        f_D[i]=0;
        //Random force
        f_R[i]=0;
        for (j=0;j<3;j++){
          f_D[i]+=P_times_dist_sqr[i][j]*(p1->m.v[j] - p2->m.v[j]);
          f_R[i]+=P_times_dist_sqr[i][j]*noise_vec[j];
        }
        //NOTE: velocity are scaled with time_step
        f_D[i]*=ia_params->dpd_pref3*omega2;
        //NOTE: noise force scales with 1/sqrt(time_step
        f_R[i]*=ia_params->dpd_pref4*omega*dist_inv;
      }
      for(j=0; j<3; j++) {
        tmp=f_R[j]-f_D[j];
        p1->f.f[j] += tmp;
        p2->f.f[j] -= tmp;
      }
  }
}

/* TODO: This function is not used anywhere. To be removed?  */
inline double inter_dpd_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
   return 0;
}
#endif

#endif

