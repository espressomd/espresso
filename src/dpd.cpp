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
#endif
