/*
  Copyright (C) 2010,2011 The ESPResSo project
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
/** \file interaction_data.c
    Implementation of interaction_data.h
 */
#include <string.h>
#include <stdlib.h>
#include "utils.h"
#include "rattle.h"
#include "interaction_data.h"
#include "errorhandling.h"
#include "communication.h"
#include "grid.h"
#include "pressure.h"
#include "p3m.h"
#include "ewald.h"
#include "debye_hueckel.h"
#include "reaction_field.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "maggs.h"
#include "elc.h"
#include "lj.h"
#include "ljgen.h"
#include "ljangle.h"
#include "steppot.h"
#include "hertzian.h"
#include "buckingham.h"
#include "soft_sphere.h"
#include "tab.h"
#include "overlap.h"
#include "ljcos.h"
#include "ljcos2.h"
#include "gb.h"
#include "parser.h"
#include "cells.h"
#include "comforce.h"
#include "comfixed.h"
#include "morse.h"
#include "dpd.h"
#include "tunable_slip.h"
#include "magnetic_non_p3m_methods.h"
#include "mdlc_correction.h"
#include "tcl_interface/interaction_data_tcl.c"

/****************************************
 * variables
 *****************************************/
int n_particle_types = 0;
int n_interaction_types = 0;
IA_parameters *ia_params = NULL;

#ifdef ADRESS
/* #ifdef THERMODYNAMIC_FORCE */
TF_parameters *tf_params = NULL;
/* #endif */
#endif

#if defined(ELECTROSTATICS) || defined(DIPOLES)
Coulomb_parameters coulomb = { 
#ifdef ELECTROSTATICS
  0.0, 0.0, COULOMB_NONE,
#endif
#ifdef DIPOLES
  0.0, 0.0, DIPOLAR_NONE,
#endif
};
#endif

#ifdef ELECTROSTATICS
Debye_hueckel_params dh_params = { 0.0, 0.0 };
Reaction_field_params rf_params = { 0.0, 0.0 };
#endif

int n_bonded_ia = 0;
Bonded_ia_parameters *bonded_ia_params = NULL;

double max_cut;
double max_cut_non_bonded;

double lj_force_cap = 0.0;
double ljangle_force_cap = 0.0;
double morse_force_cap = 0.0;
double tab_force_cap = 0.0;
double buck_force_cap = 0.0;

#ifdef CONSTRAINTS
int n_constraints       = 0;
Constraint *constraints = NULL;
#endif

/** Array containing all tabulated forces*/
DoubleList tabulated_forces;
/** Corresponding array containing all tabulated energies*/
DoubleList tabulated_energies;

#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
/** Array containing all adress tabulated forces*/
DoubleList adress_tab_forces;
/** Corresponding array containing all adress tabulated energies*/
DoubleList adress_tab_energies;
#endif

/* #ifdef THERMODYNAMIC_FORCE */
/** Array containing the thermodynamic forces **/
DoubleList thermodynamic_forces;
DoubleList thermodynamic_f_energies;
/* #endif */
#endif


/*****************************************
 * general lowlevel functions
 *****************************************/

/** Initialize force and energy tables */
void force_and_energy_tables_init() {
  init_doublelist(&tabulated_forces);
  init_doublelist(&tabulated_energies);
}

#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
/** Initialize adress force and energy tables */
void adress_force_and_energy_tables_init() {
  init_doublelist(&adress_tab_forces);
  init_doublelist(&adress_tab_energies);
}
#endif

/* #ifdef THERMODYNAMIC_FORCE */
void tf_tables_init() {
  init_doublelist(&thermodynamic_forces);
  init_doublelist(&thermodynamic_f_energies);
}
/* #endif */
#endif



/** Initialize interaction parameters. */
void initialize_ia_params(IA_parameters *params) {
#ifdef LENNARD_JONES
	params->LJ_eps =
		params->LJ_sig =
		params->LJ_cut =
		params->LJ_shift =
		params->LJ_offset =
		params->LJ_capradius = 0;
	params->LJ_min = 0;  
#endif

#ifdef LENNARD_JONES_GENERIC
  params->LJGEN_eps =
    params->LJGEN_sig =
    params->LJGEN_cut =
    params->LJGEN_shift =
    params->LJGEN_offset =
    params->LJGEN_capradius =
    params->LJGEN_a1 =
    params->LJGEN_a2 = 
    params->LJGEN_b1 =
    params->LJGEN_b2 = 0;
#endif

#ifdef LJ_ANGLE
  params->LJANGLE_eps =
    params->LJANGLE_sig =
    params->LJANGLE_cut =
    params->LJANGLE_bonded1type=
    params->LJANGLE_bonded1pos = 
    params->LJANGLE_bonded1neg = 
    params->LJANGLE_bonded2pos = 
    params->LJANGLE_bonded2neg = 
    params->LJANGLE_capradius = 0;
  params->LJANGLE_z0 = 0.;
  params->LJANGLE_dz = -1.;
  params->LJANGLE_kappa = 0.;
  params->LJANGLE_epsprime = 0.;
#endif

#ifdef SMOOTH_STEP
  params->SmSt_eps =
    params->SmSt_sig =
    params->SmSt_cut =
    params->SmSt_d =
    params->SmSt_n =
    params->SmSt_k0 = 0;
#endif

#ifdef HERTZIAN
  params->Hertzian_eps =
    params->Hertzian_sig = 0;
#endif

#ifdef BMHTF_NACL
  params->BMHTF_A =
    params->BMHTF_B =
    params->BMHTF_C =
    params->BMHTF_D =
    params->BMHTF_sig =
    params->BMHTF_cut =
    params->BMHTF_computed_shift = 0;
#endif

#ifdef MORSE
  params->MORSE_eps = 
    params->MORSE_alpha =
    params->MORSE_rmin =
    params->MORSE_cut = 
    params->MORSE_rest = 
    params->MORSE_capradius = 0;
#endif

#ifdef BUCKINGHAM
    params->BUCK_A =
    params->BUCK_B =
    params->BUCK_C =
    params->BUCK_D =
    params->BUCK_cut =
    params->BUCK_discont =
    params->BUCK_shift =
    params->BUCK_capradius = 0;
    params->BUCK_F1 = 0;
    params->BUCK_F2 = 0;
#endif

#ifdef SOFT_SPHERE
  params->soft_a =
    params->soft_n =
    params->soft_cut =
    params->soft_offset = 0;
#endif

#ifdef LJCOS
  params->LJCOS_eps =
    params->LJCOS_sig =
    params->LJCOS_cut =
    params->LJCOS_offset =
    params->LJCOS_alfa =
    params->LJCOS_beta =
    params->LJCOS_rmin = 0 ;
#endif

#ifdef LJCOS2
  params->LJCOS2_eps =
    params->LJCOS2_sig =
    params->LJCOS2_cut =
    params->LJCOS2_offset =
    params->LJCOS2_w =
    params->LJCOS2_rchange = 
    params->LJCOS2_capradius = 0 ;
#endif

#ifdef GAY_BERNE
  params->GB_eps =
    params->GB_sig =
    params->GB_cut =
    params->GB_k1 =
    params->GB_k2 =
    params->GB_mu =
    params->GB_nu =
    params->GB_chi1 =
    params->GB_chi2 = 0 ;
#endif

#ifdef TABULATED
  params->TAB_npoints = 0;
  params->TAB_startindex = 0;
  params->TAB_minval = 0.0;
  params->TAB_minval2 = 0.0;
  params->TAB_maxval = 0.0;
  params->TAB_maxval2 = 0.0;
  params->TAB_stepsize = 0.0;
  strcpy(params->TAB_filename,"");
#endif

#ifdef COMFORCE
  params->COMFORCE_flag = 0;
  params->COMFORCE_dir = 0;
  params->COMFORCE_force = 0.;
	params->COMFORCE_fratio = 0.;
#endif

#ifdef COMFIXED
  params->COMFIXED_flag = 0;
#endif

#ifdef INTER_DPD
  params->dpd_gamma = 0.0;
  params->dpd_r_cut = 0.0;
  params->dpd_wf = 0;
  params->dpd_pref1 = 0.0;
  params->dpd_pref2 = 0.0;
  params->dpd_tgamma = 0.0;
  params->dpd_tr_cut = 0.0;
  params->dpd_wf = 0;
  params->dpd_pref3 = 0;
  params->dpd_pref4 = 0;
#endif

#ifdef INTER_RF
  params->rf_on = 0;
#endif

#ifdef MOL_CUT
  params->mol_cut_type = 0;
  params->mol_cut_cutoff = 0.0;
#endif

#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
  //params->ADRESS_IC_npoints = 0;
  params->ADRESS_TAB_npoints = 0;
  params->ADRESS_TAB_startindex = 0;
  params->ADRESS_TAB_minval = 0.0;
  params->ADRESS_TAB_minval2 = 0.0;
  params->ADRESS_TAB_maxval = 0.0;
  params->ADRESS_TAB_maxval2 = 0.0;
  params->ADRESS_TAB_stepsize = 0.0;
  strcpy(params->ADRESS_TAB_filename,"");
#endif
#endif

#ifdef TUNABLE_SLIP
  params->TUNABLE_SLIP_temp  = 0.0;
  params->TUNABLE_SLIP_gamma = 0.0;
  params->TUNABLE_SLIP_r_cut = 0.0;
  params->TUNABLE_SLIP_time  = 0.0;
  params->TUNABLE_SLIP_vx  = 0.0;
  params->TUNABLE_SLIP_vy  = 0.0;
  params->TUNABLE_SLIP_vz  = 0.0;
#endif
}

#ifdef ADRESS
/* #ifdef THERMODYNAMIC_FORCE */
void initialize_tf_params(TF_parameters *params){
  params->TF_TAB_npoints = 0;
  params->TF_TAB_startindex = 0;
  
  params->TF_prefactor = 0.0;
  params->TF_TAB_minval = 0.0;
  params->TF_TAB_maxval = 0.0;
  params->TF_TAB_stepsize = 0.0;
  strcpy(params->TF_TAB_filename, "");
}
/** endif */
#endif


/** Copy interaction parameters. */
void copy_ia_params(IA_parameters *dst, IA_parameters *src) {
#ifdef LENNARD_JONES
  dst->LJ_eps = src->LJ_eps;
  dst->LJ_sig = src->LJ_sig;
  dst->LJ_cut = src->LJ_cut;
  dst->LJ_shift = src->LJ_shift;
  dst->LJ_offset = src->LJ_offset;
  dst->LJ_capradius = src->LJ_capradius;
  dst->LJ_min = src->LJ_min;
#endif

#ifdef LENNARD_JONES_GENERIC
  dst->LJGEN_eps = src->LJGEN_eps;
  dst->LJGEN_sig = src->LJGEN_sig;
  dst->LJGEN_cut = src->LJGEN_cut;
  dst->LJGEN_shift = src->LJGEN_shift;
  dst->LJGEN_offset = src->LJGEN_offset;
  dst->LJGEN_capradius = src->LJGEN_capradius;
  dst->LJGEN_a1 = src->LJGEN_a1;
  dst->LJGEN_a2 = src->LJGEN_a2;
  dst->LJGEN_b1 = src->LJGEN_b1;
  dst->LJGEN_b2 = src->LJGEN_b2;
#endif

#ifdef LJ_ANGLE
  dst->LJANGLE_eps = src->LJANGLE_eps;
  dst->LJANGLE_sig = src->LJANGLE_sig;
  dst->LJANGLE_cut = src->LJANGLE_cut;
  dst->LJANGLE_bonded1type = src->LJANGLE_bonded1type;
  dst->LJANGLE_bonded1pos = src->LJANGLE_bonded1pos;
  dst->LJANGLE_bonded1neg = src->LJANGLE_bonded1neg;
  dst->LJANGLE_bonded2pos = src->LJANGLE_bonded2pos;
  dst->LJANGLE_bonded2neg = src->LJANGLE_bonded2neg;
  dst->LJANGLE_capradius = src->LJANGLE_capradius;
  dst->LJANGLE_z0 = src->LJANGLE_z0;
  dst->LJANGLE_dz = src->LJANGLE_dz;
  dst->LJANGLE_kappa = src->LJANGLE_kappa;
  dst->LJANGLE_epsprime = src->LJANGLE_epsprime;
#endif

#ifdef SMOOTH_STEP
  dst->SmSt_eps = src->SmSt_eps;
  dst->SmSt_sig = src->SmSt_sig;
  dst->SmSt_cut = src->SmSt_cut;
  dst->SmSt_d = src->SmSt_d;
  dst->SmSt_n = src->SmSt_n;
  dst->SmSt_k0 = src->SmSt_k0;
#endif

#ifdef HERTZIAN
  dst->Hertzian_eps = src->Hertzian_eps;
  dst->Hertzian_sig = src->Hertzian_sig;
#endif

#ifdef BMHTF_NACL
  dst->BMHTF_A = src->BMHTF_A;
  dst->BMHTF_B = src->BMHTF_B;
  dst->BMHTF_C = src->BMHTF_C;
  dst->BMHTF_D = src->BMHTF_D;
  dst->BMHTF_sig = src->BMHTF_sig;
  dst->BMHTF_cut = src->BMHTF_cut;
  dst->BMHTF_computed_shift = src->BMHTF_computed_shift;
#endif

#ifdef MORSE
  dst->MORSE_eps = src->MORSE_eps;
  dst->MORSE_alpha = src->MORSE_alpha;
  dst->MORSE_rmin = src->MORSE_rmin;
  dst->MORSE_cut = src->MORSE_cut;
  dst->MORSE_rest = src->MORSE_rest;
  dst->MORSE_capradius = src->MORSE_capradius;
#endif

#ifdef BUCKINGHAM
  dst->BUCK_A = src->BUCK_A;
  dst->BUCK_B = src->BUCK_B;
  dst->BUCK_C = src->BUCK_C;
  dst->BUCK_D = src->BUCK_D;
  dst->BUCK_cut = src->BUCK_cut;
  dst->BUCK_discont = src->BUCK_discont;
  dst->BUCK_shift  = src->BUCK_shift;
  dst->BUCK_capradius = src->BUCK_capradius;
  dst->BUCK_F1 = src->BUCK_F1;
  dst->BUCK_F2 = src->BUCK_F2;
#endif

#ifdef SOFT_SPHERE
  dst->soft_a = src->soft_a;
  dst->soft_n = src->soft_n;
  dst->soft_cut = src->soft_cut;
  dst->soft_offset = src->soft_offset;
#endif

#ifdef LJCOS
  dst->LJCOS_eps = src->LJCOS_eps;
  dst->LJCOS_sig = src->LJCOS_sig;
  dst->LJCOS_cut = src->LJCOS_cut;
  dst->LJCOS_offset = src->LJCOS_offset;
  dst->LJCOS_alfa = src->LJCOS_alfa;
  dst->LJCOS_beta = src->LJCOS_beta;
  dst->LJCOS_rmin = src->LJCOS_rmin;
#endif

#ifdef LJCOS2
  dst->LJCOS2_eps       = src->LJCOS2_eps;
  dst->LJCOS2_sig       = src->LJCOS2_sig;
  dst->LJCOS2_cut       = src->LJCOS2_cut;
  dst->LJCOS2_offset    = src->LJCOS2_offset;
  dst->LJCOS2_w         = src->LJCOS2_w;
  dst->LJCOS2_rchange   = src->LJCOS2_rchange;
  dst->LJCOS2_capradius = src->LJCOS2_capradius;
#endif
  
#ifdef GAY_BERNE
  dst->GB_eps = src->GB_eps;
  dst->GB_sig = src->GB_sig;
  dst->GB_cut = src->GB_cut;
  dst->GB_k1 = src->GB_k1;
  dst->GB_k2 = src->GB_k2;
  dst->GB_mu = src->GB_mu;
  dst->GB_nu = src->GB_nu;
  dst->GB_chi1 = src->GB_chi1;
  dst->GB_chi2 = src->GB_chi2; 
#endif

#ifdef TABULATED
  dst->TAB_npoints = src->TAB_npoints;
  dst->TAB_startindex = src->TAB_startindex;
  dst->TAB_minval = src->TAB_minval;
  dst->TAB_minval2 = src->TAB_minval2;
  dst->TAB_maxval = src->TAB_maxval;
  dst->TAB_maxval2 = src->TAB_maxval2;
  dst->TAB_stepsize = src->TAB_stepsize;
  strcpy(dst->TAB_filename,src->TAB_filename);
#endif

#ifdef COMFORCE
  dst->COMFORCE_flag = src->COMFORCE_flag;
  dst->COMFORCE_dir = src->COMFORCE_dir;
  dst->COMFORCE_force = src->COMFORCE_force;
  dst->COMFORCE_fratio = src->COMFORCE_fratio;
#endif

#ifdef COMFIXED
  dst->COMFIXED_flag = src->COMFIXED_flag;
#endif

#ifdef INTER_DPD
  dst->dpd_gamma  = src->dpd_gamma;
  dst->dpd_r_cut  = src-> dpd_r_cut;
  dst->dpd_wf     = src->dpd_wf;
  dst->dpd_pref1  = src->dpd_pref1;
  dst->dpd_pref2  = src->dpd_pref2;
  dst->dpd_tgamma = src->dpd_tgamma;
  dst->dpd_tr_cut = src-> dpd_tr_cut;
  dst->dpd_twf    = src->dpd_twf;
  dst->dpd_pref3  = src->dpd_pref3;
  dst->dpd_pref4  = src->dpd_pref4;
#endif

#ifdef INTER_RF
  dst->rf_on = src->rf_on;
#endif

#ifdef MOL_CUT
  dst->mol_cut_type = src->mol_cut_type;
  dst->mol_cut_cutoff = src->mol_cut_cutoff;
#endif

#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
  //dst->ADRESS_IC_npoints = src->ADRESS_IC_npoints;
  dst->ADRESS_TAB_npoints = src->ADRESS_TAB_npoints;
  dst->ADRESS_TAB_startindex = src->ADRESS_TAB_startindex;
  dst->ADRESS_TAB_minval = src->ADRESS_TAB_minval;
  dst->ADRESS_TAB_minval2 = src->ADRESS_TAB_minval2;
  dst->ADRESS_TAB_maxval = src->ADRESS_TAB_maxval;
  dst->ADRESS_TAB_maxval2 = src->ADRESS_TAB_maxval2;
  dst->ADRESS_TAB_stepsize = src->ADRESS_TAB_stepsize;
  strcpy(dst->ADRESS_TAB_filename,src->ADRESS_TAB_filename);
#endif
#endif

#ifdef TUNABLE_SLIP
  dst->TUNABLE_SLIP_temp  = src->TUNABLE_SLIP_temp;
  dst->TUNABLE_SLIP_gamma  = src->TUNABLE_SLIP_gamma;
  dst->TUNABLE_SLIP_r_cut  = src->TUNABLE_SLIP_r_cut;
  dst->TUNABLE_SLIP_time  = src->TUNABLE_SLIP_time;
  dst->TUNABLE_SLIP_vx  = src->TUNABLE_SLIP_vx;
  dst->TUNABLE_SLIP_vy  = src->TUNABLE_SLIP_vy;
  dst->TUNABLE_SLIP_vz  = src->TUNABLE_SLIP_vz;
#endif

}

#ifdef ADRESS
/* #ifdef THERMODYNAMIC_FORCE */
void copy_tf_params(TF_parameters *dst, TF_parameters *src){
  dst->TF_TAB_npoints = src->TF_TAB_npoints;
  dst->TF_TAB_startindex = src->TF_TAB_startindex;
  dst->TF_prefactor = src->TF_prefactor;
  dst->TF_TAB_minval = src->TF_TAB_minval;
  dst->TF_TAB_maxval = src->TF_TAB_maxval;
  dst->TF_TAB_stepsize = src->TF_TAB_stepsize;
  strcpy(dst->TF_TAB_filename,src->TF_TAB_filename);
}
/* #endif */
#endif

/** returns non-zero if there is a nonbonded interaction defined */
int checkIfInteraction(IA_parameters *data) {

#ifdef LENNARD_JONES
  if (data->LJ_cut != 0)
    return 1;
#endif

#ifdef LENNARD_JONES_GENERIC
  if (data->LJGEN_cut != 0)
    return 1;
#endif

#ifdef LJ_ANGLE
  if (data->LJANGLE_cut != 0)
    return 1;
#endif

#ifdef SMOOTH_STEP
  if (data->SmSt_cut != 0)
    return 1;
#endif
  
#ifdef HERTZIAN
  if (data->Hertzian_sig != 0)
    return 1;
#endif
  
#ifdef BMHTF_NACL
  if (data->BMHTF_cut != 0)
    return 1;
#endif
  
#ifdef MORSE
  if (data->MORSE_cut != 0)
    return 1;
#endif

#ifdef BUCKINGHAM
  if (data->BUCK_cut != 0)
    return 1;
#endif

#ifdef SOFT_SPHERE
  if (data->soft_cut != 0)
    return 1;
#endif  

#ifdef LJCOS
  if (data->LJCOS_cut != 0)
    return 1;
#endif

#ifdef LJCOS2
  if (data->LJCOS2_cut != 0)
    return 1;
#endif

#ifdef GAY_BERNE
  if (data->GB_cut != 0)
    return 1;
#endif

#ifdef TABULATED
  if (data->TAB_maxval != 0)
    return 1;
#endif

#ifdef DPD
	 if (dpd_r_cut !=0)
	   return 1;
#endif

#ifdef TRANS_DPD
	 if (dpd_tr_cut !=0)
	   return 1;
#endif

#ifdef INTER_DPD
  if ( (data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0) )
    return 1;
#endif

#ifdef INTER_RF
  if (data->rf_on == 1)
    return 1;
#endif

#ifdef MOL_CUT
  if (data->mol_cut_type != 0)
    return 1;
#endif

#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
  if(data->ADRESS_TAB_maxval != 0)
    return 1;
#endif
#endif

#ifdef TUNABLE_SLIP
if (data->TUNABLE_SLIP_r_cut != 0)
    return 1;
#endif

  return 0;
}

#ifdef ADRESS
/* #ifdef THERMODYNAMIC_FORCE */
int checkIfTF(TF_parameters *data){
  if (data->TF_TAB_maxval !=0)
    return 1;
  return 0;
}
/* #endif */
#endif

char *get_name_of_bonded_ia(int i) {
  switch (i) {
  case BONDED_IA_FENE:
    return "FENE";
  case BONDED_IA_ANGLE:
    return "angle";
  case BONDED_IA_ANGLEDIST:
    return "angledist";
  case BONDED_IA_DIHEDRAL:
    return "dihedral";
  case BONDED_IA_ENDANGLEDIST:
    return "endangledist";
  case BONDED_IA_HARMONIC:
    return "HARMONIC";
  case BONDED_IA_SUBT_LJ:
    return "SUBT_LJ";
  case BONDED_IA_TABULATED:
    return "tabulated";
  case BONDED_IA_OVERLAPPED:
    return "overlapped";
  case BONDED_IA_RIGID_BOND:
    return "RIGID_BOND";
  case BONDED_IA_VIRTUAL_BOND:
    return "VIRTUAL_BOND";
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: name of unknown interaction %d requested\n",
	    this_node, i);
    errexit();
  }
  /* just to keep the compiler happy */
  return "";
}

/** This function increases the LOCAL ia_params field for non-bonded interactions
    to the given size. This function is not exported
    since it does not do this on all nodes. Use
    make_particle_type_exist for that.
*/
void realloc_ia_params(int nsize)
{
  int i, j;
  IA_parameters *new_params;
  if (nsize <= n_particle_types)
    return;

  new_params = (IA_parameters *) malloc(nsize*nsize*sizeof(IA_parameters));
  if (ia_params) {
    /* if there is an old field, copy entries and delete */
    for (i = 0; i < nsize; i++)
      for (j = 0; j < nsize; j++) {
	if ((i < n_particle_types) && (j < n_particle_types))
	  copy_ia_params(&new_params[i*nsize + j],
			 &ia_params[i*n_particle_types + j]);
	else
	  initialize_ia_params(&new_params[i*nsize + j]);
      }
    free(ia_params);
  }
  else {
    /* new field, just init */
    for (i = 0; i < nsize; i++)
      for (j = 0; j < nsize; j++)
	initialize_ia_params(&new_params[i*nsize + j]);
  }

  n_particle_types = nsize;
  ia_params = new_params;
}

#ifdef ADRESS
/* #ifdef THERMODYNAMIC_FORCE */
void realloc_tf_params(int nsize)
{
  int i;
  TF_parameters *new_params;

  if (nsize <= n_particle_types)
    return;

  new_params = (TF_parameters *) malloc(nsize*sizeof(TF_parameters));
  if (tf_params) {
    /* if there is an old field, copy entries and delete */
    for (i = 0; i < nsize; i++)
      {
	if (i < n_particle_types)
	  copy_tf_params(&new_params[i],
			 &tf_params[i]);
	else
	  initialize_tf_params(&new_params[i]);
      }
    free(tf_params);
  }
  else {
    /* new field, just init */
    for (i = 0; i < nsize; i++)
      initialize_tf_params(&new_params[i]);
  }
  
  tf_params = new_params;
}
/* #endif */
#endif

void make_particle_type_exist(int type)
{
  int ns = type + 1;
  if (ns <= n_particle_types)
    return;
  
  mpi_bcast_n_particle_types(ns);
}

void make_bond_type_exist(int type)
{
  int i, ns = type + 1;
  
  if(ns <= n_bonded_ia) {
#ifdef TABULATED
    if ( bonded_ia_params[type].type == BONDED_IA_TABULATED && 
	 bonded_ia_params[type].p.tab.npoints > 0 ) {
      free(bonded_ia_params[type].p.tab.f);
      free(bonded_ia_params[type].p.tab.e);
    }
#endif 
#ifdef OVERLAPPED
    if ( bonded_ia_params[type].type == BONDED_IA_OVERLAPPED &&
         bonded_ia_params[type].p.overlap.noverlaps > 0 ) {
      free(bonded_ia_params[type].p.overlap.para_a);
      free(bonded_ia_params[type].p.overlap.para_b);
      free(bonded_ia_params[type].p.overlap.para_c);
    }
#endif
    return;
  }
  /* else allocate new memory */
  bonded_ia_params = (Bonded_ia_parameters *)realloc(bonded_ia_params,
						     ns*sizeof(Bonded_ia_parameters));
  /* set bond types not used as undefined */
  for (i = n_bonded_ia; i < ns; i++)
    bonded_ia_params[i].type = BONDED_IA_NONE;
 
  n_bonded_ia = ns;
}

void calc_maximal_cutoff()
{
  int i, j;
  double max_cut_tmp;
  double max_cut_bonded=-1.0;
  max_cut = -1.0;
  max_cut_non_bonded = -1.0;

  /* bonded */
  for (i = 0; i < n_bonded_ia; i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_FENE:
      max_cut_tmp = bonded_ia_params[i].p.fene.r0+bonded_ia_params[i].p.fene.drmax;
      if(max_cut_bonded < max_cut_tmp)
	max_cut_bonded = max_cut_tmp;
      break;
    case BONDED_IA_HARMONIC:
      if((bonded_ia_params[i].p.harmonic.r_cut>0)&&(max_cut_bonded < bonded_ia_params[i].p.harmonic.r_cut))
	max_cut_bonded = bonded_ia_params[i].p.harmonic.r_cut;
      break;
    case BONDED_IA_SUBT_LJ:
      if(max_cut_bonded < bonded_ia_params[i].p.subt_lj.r)
	max_cut_bonded = bonded_ia_params[i].p.subt_lj.r;
      break;
    case BONDED_IA_RIGID_BOND:
      if(max_cut_bonded < sqrt(bonded_ia_params[i].p.rigid_bond.d2))
	max_cut_bonded = sqrt(bonded_ia_params[i].p.rigid_bond.d2);
       break;
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      if(bonded_ia_params[i].p.tab.type == TAB_BOND_LENGTH &&
	 max_cut_bonded < bonded_ia_params[i].p.tab.maxval)
	max_cut_bonded = bonded_ia_params[i].p.tab.maxval;
      break;
#endif
#ifdef OVERLAPPED 
    case BONDED_IA_OVERLAPPED:
      /* in UNIT Angstrom */
      if(bonded_ia_params[i].p.overlap.type == OVERLAP_BOND_LENGTH &&
         max_cut_bonded < bonded_ia_params[i].p.overlap.maxval)
        max_cut_bonded = bonded_ia_params[i].p.overlap.maxval;
      break;
#endif
    default:
     break;
    }
  }

  /* Bond angle and dihedral potentials do not contain a cutoff
     intrinsically. The cutoff for these potentials depends on the
     bond length potentials. For bond angle potentials nothing has to
     be done (it is assumed, that particles participating in a bond
     angle or dihedral potential are bound to each other by some bond
     length potential (FENE, Harmonic or tabulated)). For dihedral
     potentials (both normal and tabulated ones) it follows, that the
     cutoff is TWO TIMES the maximal cutoff! That's what the following
     lines assure. */
  max_cut_tmp = 2.0*max_cut_bonded;
  for (i = 0; i < n_bonded_ia; i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_DIHEDRAL:
      max_cut_bonded = max_cut_tmp;
      break;
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      if(bonded_ia_params[i].p.tab.type == TAB_BOND_DIHEDRAL)
	max_cut_bonded = max_cut_tmp;
      break;
#endif
#ifdef OVERLAPPED 
    case BONDED_IA_OVERLAPPED:
      if(bonded_ia_params[i].p.overlap.type == OVERLAP_BOND_DIHEDRAL)
        max_cut_bonded = max_cut_tmp;
      break;
#endif
    default:
      break;
    }
  }
  max_cut=max_cut_bonded;

  /* non bonded */
  for (i = 0; i < n_particle_types; i++)
     for (j = i; j < n_particle_types; j++) {
       if (checkIfParticlesInteract(i, j)) {
	 IA_parameters *data = get_ia_param(i, j);
#ifdef LENNARD_JONES
	 if (data->LJ_cut != 0) {
	   if(max_cut_non_bonded < (data->LJ_cut+data->LJ_offset) )
	     max_cut_non_bonded = (data->LJ_cut+data->LJ_offset);
	 }
#endif

#ifdef DPD
	 if (dpd_r_cut !=0) {
	   if(max_cut_non_bonded < dpd_r_cut)
	     max_cut_non_bonded = dpd_r_cut;
	 }
#endif

#ifdef TRANS_DPD
	 if (dpd_tr_cut !=0) {
	   if(max_cut_non_bonded < dpd_tr_cut)
	     max_cut_non_bonded = dpd_tr_cut;
	 }
#endif

#ifdef LENNARD_JONES_GENERIC
	 if (data->LJGEN_cut != 0) {
	   if(max_cut_non_bonded < (data->LJGEN_cut+data->LJGEN_offset) )
	     max_cut_non_bonded = (data->LJGEN_cut+data->LJGEN_offset);
	 }
#endif

#ifdef LJ_ANGLE
	 if (data->LJANGLE_cut != 0) {
	   if(max_cut_non_bonded < (data->LJANGLE_cut) )
	     max_cut_non_bonded = (data->LJANGLE_cut);
	 }
#endif

#ifdef INTER_DPD
	 if ((data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0)){
	   if(max_cut_non_bonded < ( (data->dpd_r_cut > data->dpd_tr_cut)?data->dpd_r_cut:data->dpd_tr_cut ) )
	     max_cut_non_bonded = ( (data->dpd_r_cut > data->dpd_tr_cut)?data->dpd_r_cut:data->dpd_tr_cut );
	 }
#endif

#ifdef SMOOTH_STEP
         if (data->SmSt_cut != 0) {
           if(max_cut_non_bonded < data->SmSt_cut)
             max_cut_non_bonded = data->SmSt_cut;
         }
#endif

#ifdef HERTZIAN
         if (data->Hertzian_sig != 0) {
           if(max_cut_non_bonded < data->Hertzian_sig)
             max_cut_non_bonded = data->Hertzian_sig;
         }
#endif

#ifdef BMHTF_NACL
         if (data->BMHTF_cut != 0) {
           if(max_cut_non_bonded < data->BMHTF_cut)
             max_cut_non_bonded = data->BMHTF_cut;
         }
#endif

#ifdef MORSE
         if (data->MORSE_cut != 0) {
           if(max_cut_non_bonded < (data->MORSE_cut) )
             max_cut_non_bonded = (data->MORSE_cut);
         }
#endif

#ifdef BUCKINGHAM
	 if (data->BUCK_cut != 0) {
	   if(max_cut_non_bonded < data->BUCK_cut )
	     max_cut_non_bonded = data->BUCK_cut;
	 }
#endif

#ifdef SOFT_SPHERE
	 if (data->soft_cut != 0) {
	   if(max_cut_non_bonded < data->soft_cut )
	     max_cut_non_bonded = data->soft_cut;
	 }
#endif

#ifdef LJCOS
	 if (data->LJCOS_cut != 0) {
	   if(max_cut_non_bonded < (data->LJCOS_cut+data->LJCOS_offset) )
	     max_cut_non_bonded = (data->LJCOS_cut+data->LJCOS_offset);
	 }
#endif

#ifdef LJCOS2
	 if (data->LJCOS2_cut != 0) {
	   if(max_cut_non_bonded < (data->LJCOS2_cut+data->LJCOS2_offset) )
	     max_cut_non_bonded = (data->LJCOS2_cut+data->LJCOS2_offset);
	 }
#endif

#ifdef GAY_BERNE
	 if (data->GB_cut != 0) {
	   if(max_cut_non_bonded < (data->GB_cut) )
	     max_cut_non_bonded = (data->GB_cut);
	 }
#endif

#ifdef TABULATED
	 if (data->TAB_maxval != 0){
	   if(max_cut_non_bonded < (data->TAB_maxval ))
	     max_cut_non_bonded = data->TAB_maxval;
	 }
#endif
	 
#ifdef ADRESS
#ifdef INTERFACE_CORRECTION
	 if (data->ADRESS_TAB_maxval !=0){
	   if(max_cut_non_bonded < (data->ADRESS_TAB_maxval ))
	     max_cut_non_bonded = data->ADRESS_TAB_maxval;
	 }
#endif
#endif

#ifdef TUNABLE_SLIP
	 if (data->TUNABLE_SLIP_r_cut != 0){
	   if(max_cut_non_bonded < (data->TUNABLE_SLIP_r_cut ))
	     max_cut_non_bonded = data->TUNABLE_SLIP_r_cut;
	 }
#endif
       }
     }

#ifdef ELECTROSTATICS
  /* real space electrostatic */
  switch (coulomb.method) {
#ifdef P3M 
  case COULOMB_ELC_P3M:
    if (max_cut_non_bonded < elc_params.space_layer)
      max_cut_non_bonded = elc_params.space_layer;
    // fall through
  case COULOMB_P3M:
    if (max_cut_non_bonded < p3m.params.r_cut)
      max_cut_non_bonded = p3m.params.r_cut;
    break;
#endif
  case COULOMB_EWALD:
    if (max_cut_non_bonded < ewald.r_cut)
      max_cut_non_bonded = ewald.r_cut;
    break;
  case COULOMB_DH:
    if (max_cut_non_bonded < dh_params.r_cut)
      max_cut_non_bonded = dh_params.r_cut;
    break;
  case COULOMB_RF:
  case COULOMB_INTER_RF:
    if (max_cut_non_bonded < rf_params.r_cut)
      max_cut_non_bonded = rf_params.r_cut;
    break;
  case COULOMB_MMM1D:
    /* needs n-squared calculation anyways */
    if (max_cut_non_bonded < 0)
      max_cut_non_bonded = 0;
    break;
  case COULOMB_MMM2D:
    /* needs n-squared rsp. layered calculation, and
       it is pretty complicated to find the minimal
       required cell height. */
    if (max_cut_non_bonded < 0)
      max_cut_non_bonded = 0;
    break;
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DP3M
  switch (coulomb.Dmethod) {
  case DIPOLAR_P3M:
    if (max_cut_non_bonded < dp3m.params.r_cut)
      max_cut_non_bonded = dp3m.params.r_cut;
    break;
  }       
#endif /*ifdef DP3M */

#ifdef MOL_CUT
if(max_cut_bonded > 0)
	max_cut_non_bonded +=2.0* max_cut_bonded;
  else max_cut_non_bonded=2.0* max_cut_bonded;
#ifdef ONE_PROC_ADRESS
  max_cut_non_bonded -= 2.0*max_cut_bonded;
#endif
#endif
  /* make max_cut the maximal cutoff of both bonded and non-bonded interactions */
  if ( max_cut_non_bonded > max_cut) max_cut = max_cut_non_bonded;
}

int check_obs_calc_initialized()
{
  /* set to zero if initialization was not successful. */
  int state = 1;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_MMM1D: if (MMM1D_sanity_checks()) state = 0; break;
  case COULOMB_MMM2D: if (MMM2D_sanity_checks()) state = 0; break;
#ifdef P3M
  case COULOMB_ELC_P3M: if (ELC_sanity_checks()) state = 0; // fall through
  case COULOMB_P3M: if (p3m_sanity_checks()) state = 0; break;
#endif
  case COULOMB_EWALD: if (EWALD_sanity_checks()) state = 0; break;
  }
#endif /* ifdef ELECTROSTATICS */

#ifdef DIPOLES
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M: if (mdlc_sanity_checks()) state = 0; // fall through
  case DIPOLAR_P3M: if (dp3m_sanity_checks()) state = 0; break;
#endif
  case DIPOLAR_MDLC_DS: if (mdlc_sanity_checks()) state = 0; // fall through
  case DIPOLAR_DS: if (magnetic_dipolar_direct_sum_sanity_checks()) state = 0; break;
  }
#endif /* ifdef  DIPOLES */

  return state;
}


#ifdef ELECTROSTATICS

/********************************************************************************/
/*                                 electrostatics                               */
/********************************************************************************/

int coulomb_set_bjerrum(double bjerrum)
{
  if (bjerrum < 0.0)
    return TCL_ERROR;
  
  coulomb.bjerrum = bjerrum;

  if (coulomb.bjerrum == 0.0) {
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_ELC_P3M:
    case COULOMB_P3M:
      p3m_set_bjerrum();
      break;
#endif
    case COULOMB_EWALD:
      ewald.alpha    = 0.0;
      ewald.alpha_L  = 0.0;
      ewald.r_cut    = 0.0;
      ewald.r_cut_iL = 0.0;
      break;
    case COULOMB_DH:
      dh_params.r_cut   = 0.0;
      dh_params.kappa   = 0.0;
    case COULOMB_RF:
    case COULOMB_INTER_RF:
      rf_params.kappa  = 0.0;
      rf_params.epsilon1   = 0.0;
      rf_params.epsilon2   = 0.0;
      rf_params.r_cut   = 0.0;
      rf_params.B   = 0.0;
    case COULOMB_MMM1D:
      mmm1d_params.maxPWerror = 1e40;
      mmm1d_params.bessel_cutoff = 0;
    }
 
    mpi_bcast_coulomb_params();
    coulomb.method = COULOMB_NONE;
    mpi_bcast_coulomb_params();

  }

  return TCL_OK;
}

#endif /*ifdef ELECTROSTATICS */


#ifdef DIPOLES
int dipolar_set_Dbjerrum(double bjerrum)
{
  if (bjerrum < 0.0)
    return TCL_ERROR;
  
  coulomb.Dbjerrum = bjerrum;

  if (coulomb.Dbjerrum == 0.0) {
    switch (coulomb.Dmethod) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:
    case DIPOLAR_P3M:
      coulomb.Dbjerrum = bjerrum;
      dp3m_set_bjerrum();
      break;
#endif
    }
 
    mpi_bcast_coulomb_params();
    coulomb.Dmethod = DIPOLAR_NONE;
    mpi_bcast_coulomb_params();

  }

  return TCL_OK;
}
#endif


/********************************************************************************/
/*                                       parsing                                */
/********************************************************************************/

#ifdef BOND_VIRTUAL
int virtual_set_params(int bond_type)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].type = BONDED_IA_VIRTUAL_BOND;
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}
#endif

int lennard_jones_set_params(int part_type_a, int part_type_b,
				      double eps, double sig, double cut,
				      double shift, double offset,
				      double cap_radius, double min)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* LJ should be symmetrically */
  data->LJ_eps    = data_sym->LJ_eps    = eps;
  data->LJ_sig    = data_sym->LJ_sig    = sig;
  data->LJ_cut    = data_sym->LJ_cut    = cut;
  data->LJ_shift  = data_sym->LJ_shift  = shift;
  data->LJ_offset = data_sym->LJ_offset = offset;
 
  if (cap_radius > 0) {
    data->LJ_capradius = cap_radius;
    data_sym->LJ_capradius = cap_radius;
  }

  if (min > 0) {
	  data->LJ_min = min;
	  data_sym->LJ_min = min;	  
  }
  

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (lj_force_cap != -1.0)
    mpi_lj_cap_forces(lj_force_cap);

  return TCL_OK;
}
