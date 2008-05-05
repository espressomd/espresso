// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while u tilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2007; all rights reserved unless otherwise stated.
#ifndef DPD_H
#define DPD_H
/** \file dpd.h
 *  Routines to use dpd as thermostat or pair force
 *  T. Soddemann, B. Duenweg and K. Kremer, Phys. Rev. E 68, 046702 (2003)
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:junghans@mpip-mainz.mpg.de">Christoph</a>
*/

#include <tcl.h>
#include "thermostat.h"
#include "statistics.h"

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

void dpd_parse_off(Tcl_Interp *interp, int argc, char **argv);
int thermo_parse_dpd(Tcl_Interp *interp, int argc, char **argv);
void dpd_print(Tcl_Interp *interp);
void thermo_init_dpd();
void dpd_usage(Tcl_Interp *interp, int argc, char **argv);
void dpd_heat_up();
void dpd_cool_down();

/** Calculate Random Force and Friction Force acting between particle
    p1 and p2 and add them to their forces. */
MDINLINE void add_dpd_thermo_pair_force(Particle *p1, Particle *p2, double d[3], double dist, double dist2)
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

#ifdef DPD_MASS_RED
  massf=2*PMASS(*p1)*PMASS(*p2)/(PMASS(*p1)+PMASS(*p2));
#endif

#ifdef DPD_MASS_LIN
  massf=0.5*(PMASS(*p1)+PMASS(*p2));
#endif

  double com_dist;
  //double intra_fractor;

#ifdef MOL_CUT
  com_dist=get_mol_dist(p1,p2);
#else
  com_dist=dist;
#endif

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
    //omega*=intra_fractor;
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
  if ((com_dist < dpd_tr_cut)&&(dpd_tgamma > 0.0)){
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
      //omega*=intra_fractor;
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
void interdpd_heat_up();
void interdpd_cool_down();
void interdpd_parse_off();
int printinterdpdIAToResult(Tcl_Interp *interp, int i, int j);
int interdpd_set_params(int part_type_a, int part_type_b,
				      double gamma, double r_c, int wf,
				      double tgamma, double tr_c,
				      int twf);
int thermo_parse_interdpd(Tcl_Interp *interp, int argc, char ** argv);
int interdpd_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv);
void interdpd_init();
void interdpd_update_params(double pref2_scale);

MDINLINE void add_interdpd_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
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
  double P_times_dist_sqr[3][3]={{dist2,0,0},{0,dist2,0},{0,0,dist2}},noise_vec[3];
  double f_D[3],f_R[3];
  double tmp;
#ifdef DPD_MASS
  double massf;
#endif

#ifdef DPD_MASS_RED
  massf=2*PMASS(*p1)*PMASS(*p2)/(PMASS(*p1)+PMASS(*p2));
#endif

#ifdef DPD_MASS_LIN
  massf=0.5*(PMASS(*p1)+PMASS(*p2));
#endif

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


MDINLINE double interdpd_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
   return 0;
}
#endif

#endif

