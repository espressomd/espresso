/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file dpd.cpp
    Implementation of \ref dpd.hpp "dpd.hpp"
 */
#include "dpd.hpp"

/** Flag to decide wether to allow for fixed particles with DPD */
int dpd_ignore_fixed_particles=1;

/* DPD THERMOSTAT */
/* DPD longitudinal friction coefficient gamma. */
double dpd_gamma = 0.0;
/* DPD thermostat cutoff */
double dpd_r_cut = 0.0;
/* DPD weightfunction */
int dpd_wf = 0;

/* DPD transversal friction coefficient gamma. */
double dpd_tgamma = 0.0;
/* DPD thermostat trans cutoff */
double dpd_tr_cut = 0.0;
/* trans DPD weightfunction */
int dpd_twf = 0;

#ifdef DPD

/** Chatterjee 2007 proposes that for DPD with Lees Edwards BCs,
 *  it is better not to count interactions with Ghost particles. */
static bool le_chatterjee_test_pair(Particle *p1, Particle *p2){
#ifdef LEES_EDWARDS
    /* cannot assume that we have a flag set to mark ghost particles,
     * but can assume (for LE) that non-ghost
     * particle y-coords are always imaged inside the box */

    if( p1->r.p[1] < 0 )        return true;
    if( p1->r.p[1] > box_l[1] ) return true;
    if( p2->r.p[1] < 0 )        return true;
    if( p2->r.p[1] > box_l[1] ) return true;
#endif
    return false;
}

/* inverse off DPD thermostat cutoff */
double dpd_r_cut_inv = 0.0;
double dpd_pref1;
double dpd_pref2;
static double dpd_pref2_buffer;

#ifdef TRANS_DPD 
/* inverse off trans DPD thermostat cutoff */
double dpd_tr_cut_inv = 0.0;
double dpd_pref3;
double dpd_pref4;
static double dpd_pref4_buffer;
#endif

void dpd_switch_off()
{
  extern double dpd_gamma,dpd_r_cut;
  extern int dpd_wf;
#ifdef TRANS_DPD
  extern double dpd_tgamma,dpd_tr_cut;
  extern int dpd_twf;
#endif
  dpd_gamma = 0;
  mpi_bcast_parameter(FIELD_DPD_GAMMA);
  dpd_r_cut = 0;
  mpi_bcast_parameter(FIELD_DPD_RCUT);
  dpd_wf=0;
  mpi_bcast_parameter(FIELD_DPD_WF);
#ifdef TRANS_DPD
  dpd_tgamma = 0;
  mpi_bcast_parameter(FIELD_DPD_TGAMMA);
  dpd_tr_cut=0;
  mpi_bcast_parameter(FIELD_DPD_TRCUT);
  dpd_twf=0;
  mpi_bcast_parameter(FIELD_DPD_TWF);
#endif
}


void thermo_init_dpd()
{
  extern double dpd_gamma,dpd_r_cut,dpd_pref1,dpd_pref2;
  /*extern int dpd_wf;*/
#ifdef TRANS_DPD
  extern double dpd_tgamma,dpd_tr_cut,dpd_pref3,dpd_pref4;
  /*extern int dpd_twf;*/
#endif
  /* prefactor friction force */
  /* NOTE: velocities are scaled with time_step, so divide by time_step here*/
  dpd_pref1 = dpd_gamma/time_step;  
  /* prefactor random force */
  /*NOTE random force is propto sqrt(time_step)*/
  dpd_pref2 = sqrt(24.0*temperature*dpd_gamma/time_step);
  dpd_r_cut_inv = 1.0/dpd_r_cut;
#ifdef TRANS_DPD
  /* NOTE: velocities are scaled with time_step, so divide by time_step here*/
  dpd_pref3 = dpd_tgamma/time_step;
  /*NOTE random force is propto sqrt(time_step)*/
  dpd_pref4 = sqrt(24.0*temperature*dpd_tgamma/time_step);
  dpd_tr_cut_inv = 1.0/dpd_tr_cut;
#endif
  THERMO_TRACE(fprintf(stderr,"%d: thermo_init_dpd: dpd_pref1=%f, dpd_pref2=%f",
		       this_node,dpd_pref1,dpd_pref2));
#ifdef TRANS_DPD
  THERMO_TRACE(fprintf(stderr,",dpd_pref3=%f, dpd_pref4=%f\n",dpd_pref3,dpd_pref4));
#endif
  THERMO_TRACE(fprintf(stderr,"\n"));
}

void dpd_heat_up()
{
   extern double dpd_pref2;
   extern double dpd_pref2_buffer;
#ifdef TRANS_DPD
   extern double dpd_pref4;
   extern double dpd_pref4_buffer;
#endif
      dpd_pref2_buffer = dpd_pref2;
      dpd_pref2 *= sqrt(3);
#ifdef TRANS_DPD
      dpd_pref4_buffer = dpd_pref4;
      dpd_pref4 *= sqrt(3);
#endif
}


void dpd_cool_down()
{
   extern double dpd_pref2;
   extern double dpd_pref2_buffer;
#ifdef TRNAS_DPD
   extern double dpd_pref4;
   extern double dpd_pref4_buffer;
#endif
      dpd_pref2 = dpd_pref2_buffer;
#ifdef TRANS_DPD
      dpd_pref4 = dpd_pref4_buffer;
#endif
}

void add_dpd_thermo_pair_force(Particle * p1, Particle * p2, double d[3], double dist, double dist2)
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

  if( le_chatterjee_test_pair(p1, p2) ) return;
  
#ifdef EXTERNAL_FORCES
  // if any of the two particles is fixed in some direction then
  // do not add any dissipative or stochastic dpd force part
  // because dissipation-fluctuation theorem is violated
  if (dpd_ignore_fixed_particles)
    if ( (p1->p.ext_flag | p2->p.ext_flag) & COORDS_FIX_MASK) return;
#endif

#ifdef VIRTUAL_SITES
    if (ifParticleIsVirtual(p1) || ifParticleIsVirtual(p2)) return;
#endif

#ifdef DPD_MASS_RED
  massf=2*(*p1).p.mass*(*p2).p.mass/((*p1).p.mass+(*p2).p.mass);
#endif

#ifdef DPD_MASS_LIN
  massf=0.5*((*p1).p.mass+(*p2).p.mass);
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
        //NOTE: noise force scales with 1/sqrt(time_step)
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
void inter_dpd_heat_up()
{
	double pref_scale=sqrt(3);
	inter_dpd_update_params(pref_scale);
}


void inter_dpd_cool_down()
{
	double pref_scale=1.0/sqrt(3);
	inter_dpd_update_params(pref_scale);
}

int inter_dpd_set_params(int part_type_a, int part_type_b,
			 double gamma, double r_c, int wf,
			 double tgamma, double tr_c,
			 int twf)
{
  extern double temperature;
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  data->dpd_gamma  = gamma;
  data->dpd_r_cut  = r_c;
  data->dpd_wf     = wf;
  data->dpd_pref1  = gamma/time_step;
  data->dpd_pref2  = sqrt(24.0*temperature*gamma/time_step);
  data->dpd_tgamma = tgamma;
  data->dpd_tr_cut = tr_c;
  data->dpd_twf    = twf;
  data->dpd_pref3  = tgamma/time_step;
  data->dpd_pref4  = sqrt(24.0*temperature*tgamma/time_step);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

void inter_dpd_init(){
   extern double temperature;
   int type_a,type_b;
   IA_parameters *data;

   for (type_a=0;type_a<n_particle_types;type_a++){
      for (type_b=0;type_b<n_particle_types;type_b++){
         data=get_ia_param(type_a,type_b);
         if ( (data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0) ) {
            data->dpd_pref1=data->dpd_gamma/time_step;
            data->dpd_pref2=sqrt(24.0*temperature*data->dpd_gamma/time_step);
            data->dpd_pref3=data->dpd_tgamma/time_step;
            data->dpd_pref4=sqrt(24.0*temperature*data->dpd_tgamma/time_step);
         }
      }
   }
}

void inter_dpd_switch_off(void){
   int type_a,type_b;
   IA_parameters *data;
   for (type_a=0;type_a<n_particle_types;type_a++){
      for (type_b=0;type_b<n_particle_types;type_b++){
         data=get_ia_param(type_a,type_b);
         data->dpd_gamma  = data->dpd_r_cut  = data->dpd_wf =
         data->dpd_pref1  = data->dpd_pref2  = data->dpd_tgamma =
         data->dpd_tr_cut = data->dpd_twf    = data->dpd_pref3  =
         data->dpd_pref4  = 0.0;
      }
   }
}

void inter_dpd_update_params(double pref_scale)
{
   int type_a,type_b;
   IA_parameters *data;

   for (type_a=0;type_a<n_particle_types;type_a++){
      for (type_b=0;type_b<n_particle_types;type_b++){
         data=get_ia_param(type_a,type_b);
         if ( (data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0) ) {
            data->dpd_pref2*=pref_scale;
            data->dpd_pref4*=pref_scale;
         }
      }
   }
}

void add_inter_dpd_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
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
  // Prohibits calculation of velocity dependent
  // force between two fixed particles.
  if (p1->p.ext_flag & p2->p.ext_flag & COORDS_FIX_MASK) 
    return;
  // if any of the two particles is fixed in some direction then
  // do not add any dissipative or stochastic dpd force part
  // because dissipation-fluctuation theorem is violated
  if (dpd_ignore_fixed_particles)
    if ( (p1->p.ext_flag | p2->p.ext_flag) & COORDS_FIX_MASK) return;
#endif

#ifdef DPD
  if( le_chatterjee_test_pair(p1, p2) ) return;
#endif  
#ifdef DPD_MASS_RED
  massf=2*(*p1).p.mass*(*p2).p.mass/((*p1).p.mass+(*p2).p.mass);
#endif

#ifdef DPD_MASS_LIN
  massf=0.5*((*p1).p.mass+(*p2).p.mass);
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

#endif
