// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
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
#include "steppot.h"
#include "buckingham.h"
#include "soft_sphere.h"
#include "tab.h"
#include "ljcos.h"
#include "ljcos2.h"
#include "gb.h"
#include "parser.h"
#include "cells.h"
#include "comforce.h"
#include "comfixed.h"
#include "morse.h"
#include "dpd.h"

/****************************************
 * variables
 *****************************************/
int n_particle_types = 0;
int n_interaction_types = 0;
IA_parameters *ia_params = NULL;

#ifdef ELECTROSTATICS
Coulomb_parameters coulomb = { 0.0, 0.0, COULOMB_NONE };
Debye_hueckel_params dh_params = { 0.0, 0.0 };
Reaction_field_params rf_params = { 0.0, 0.0 };
#endif

int n_bonded_ia = 0;
Bonded_ia_parameters *bonded_ia_params = NULL;

double max_cut;
double max_cut_non_bonded;

double lj_force_cap = 0.0;
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

///
int printCoulombIAToResult(Tcl_Interp *interp);

/*****************************************
 * general lowlevel functions
 *****************************************/

/** Initialize force and energy tables */
void force_and_energy_tables_init() {
  init_doublelist(&tabulated_forces);
  init_doublelist(&tabulated_energies);
}

/** Initialize interaction parameters. */
void initialize_ia_params(IA_parameters *params) {
#ifdef LENNARD_JONES
  params->LJ_eps =
    params->LJ_sig =
    params->LJ_cut =
    params->LJ_shift =
    params->LJ_offset =
    params->LJ_capradius = 0;
#endif

#ifdef LENNARD_JONES_GENERIC
  params->LJGEN_eps =
    params->LJGEN_sig =
    params->LJGEN_cut =
    params->LJGEN_shift =
    params->LJGEN_offset =
    params->LJGEN_capradius =
    params->LJGEN_a1 =
    params->LJGEN_a2 = 0;
#endif

#ifdef SMOOTH_STEP
  params->SmSt_eps =
    params->SmSt_sig =
    params->SmSt_cut =
    params->SmSt_d =
    params->SmSt_n =
    params->SmSt_k0 = 0;
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

#ifdef ROTATION
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
  params->dpd_temp = -1;
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
  params->rf_coul_pref = 0;
  params->rf_kappa = 0;
  params->rf_epsilon1 = 0;
  params->rf_epsilon2 = 0;
  params->rf_r_cut = 0;
  params->rf_B0 = 0;
  params->rf_B1 = 0;
#endif
}

/** Copy interaction parameters. */
void copy_ia_params(IA_parameters *dst, IA_parameters *src) {
#ifdef LENNARD_JONES
  dst->LJ_eps = src->LJ_eps;
  dst->LJ_sig = src->LJ_sig;
  dst->LJ_cut = src->LJ_cut;
  dst->LJ_shift = src->LJ_shift;
  dst->LJ_offset = src->LJ_offset;
  dst->LJ_capradius = src->LJ_capradius;
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
#endif

#ifdef SMOOTH_STEP
  dst->SmSt_eps = src->SmSt_eps;
  dst->SmSt_sig = src->SmSt_sig;
  dst->SmSt_cut = src->SmSt_cut;
  dst->SmSt_d = src->SmSt_d;
  dst->SmSt_n = src->SmSt_n;
  dst->SmSt_k0 = src->SmSt_k0;
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
  
#ifdef ROTATION
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
  dst->dpd_temp   = src->dpd_temp;
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
  dst->rf_coul_pref = src->rf_coul_pref;
  dst->rf_kappa     = src->rf_kappa;
  dst->rf_epsilon1  = src->rf_epsilon1;
  dst->rf_epsilon2  = src->rf_epsilon2;
  dst->rf_r_cut     = src->rf_r_cut;
  dst->rf_B0        = src->rf_B0;
  dst->rf_B1        = src->rf_B1;
#endif

}

/** returns non-zero if particles of type i and j have a nonbonded interaction */
int checkIfParticlesInteract(int i, int j) {
  IA_parameters *data = get_ia_param(i, j);

#ifdef LENNARD_JONES
  if (data->LJ_cut != 0)
    return 1;
#endif

#ifdef LENNARD_JONES_GENERIC
  if (data->LJGEN_cut != 0)
    return 1;
#endif

#ifdef SMOOTH_STEP
  if (data->SmSt_cut != 0)
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

#ifdef ROTATION
  if (data->GB_cut != 0)
    return 1;
#endif

#ifdef TABULATED
  if (data->TAB_maxval != 0)
    return 1;
#endif

#ifdef INTER_DPD
  if ( (data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0) )
    return 1;
#endif

#ifdef INTER_RF
  if (data->rf_r_cut != 0)
    return 1;
#endif

  return 0;
}

char *get_name_of_bonded_ia(int i) {
  switch (i) {
  case BONDED_IA_FENE:
    return "FENE";
  case BONDED_IA_ANGLE:
    return "angle";
  case BONDED_IA_DIHEDRAL:
    return "dihedral";
  case BONDED_IA_HARMONIC:
    return "HARMONIC";
  case BONDED_IA_SUBT_LJ:
    return "SUBT_LJ";
  case BONDED_IA_TABULATED:
    return "tabulated";
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
  max_cut = -1.0;
  max_cut_non_bonded = -1.0;

  /* bonded */
  for (i = 0; i < n_bonded_ia; i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_FENE:
      if(max_cut < bonded_ia_params[i].p.fene.r)
	max_cut = bonded_ia_params[i].p.fene.r;
      break;
    case BONDED_IA_HARMONIC:
      if(max_cut < bonded_ia_params[i].p.harmonic.r)
	max_cut = bonded_ia_params[i].p.harmonic.r;
      break;
    case BONDED_IA_SUBT_LJ:
      if(max_cut < bonded_ia_params[i].p.subt_lj.r)
	max_cut = bonded_ia_params[i].p.subt_lj.r;
      break;
    case BONDED_IA_RIGID_BOND:
      if(max_cut < sqrt(bonded_ia_params[i].p.rigid_bond.d2))
	max_cut = sqrt(bonded_ia_params[i].p.rigid_bond.d2);
       break;
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      if(bonded_ia_params[i].p.tab.type == TAB_BOND_LENGTH &&
	 max_cut < bonded_ia_params[i].p.tab.maxval)
	max_cut = bonded_ia_params[i].p.tab.maxval;
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
  max_cut_tmp = 2.0*max_cut;
  for (i = 0; i < n_bonded_ia; i++) {
    switch (bonded_ia_params[i].type) {
    case BONDED_IA_DIHEDRAL:
      max_cut = max_cut_tmp;
      break;
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      if(bonded_ia_params[i].p.tab.type == TAB_BOND_DIHEDRAL)
	max_cut = max_cut_tmp;
      break;
#endif
    default:
      break;
    }
  }

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

#ifdef INTER_RF
	 if (data->rf_r_cut != 0) {
	   if(max_cut_non_bonded < data->rf_r_cut )
	     max_cut_non_bonded = data->rf_r_cut;
	 }
#endif

#ifdef ROTATION
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
       }
     }
#ifdef ELECTROSTATICS
  /* real space electrostatic */
  switch (coulomb.method) {
#ifdef ELP3M
  case COULOMB_ELC_P3M:
    if (max_cut_non_bonded < elc_params.space_layer)
      max_cut_non_bonded = elc_params.space_layer;
    // fall through
  case COULOMB_P3M:
    if (max_cut_non_bonded < p3m.r_cut)
      max_cut_non_bonded = p3m.r_cut;
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
  case COULOMB_MAGGS:
    /* for the Yukawa subtraction scheme 
       the real cut off is needed */
    if((maggs.yukawa == 1) && (max_cut_non_bonded < maggs.r_cut))
      max_cut_non_bonded = maggs.r_cut;
    /* the cell size depends also on the grid spacing */
    if((maggs.yukawa == 0) && (max_cut < maggs.a))
      max_cut = maggs.a;  
    break;
  }
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
#ifdef ELP3M
  case COULOMB_ELC_P3M: if (ELC_sanity_checks()) state = 0; // fall through
  case COULOMB_P3M: if (P3M_sanity_checks()) state = 0; break;
#endif
  case COULOMB_EWALD: if (EWALD_sanity_checks()) state = 0; break;
  case COULOMB_MAGGS: if (Maggs_sanity_checks()) state = 0; break;
  }
#endif

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
#ifdef ELP3M
    case COULOMB_ELC_P3M:
    case COULOMB_P3M:
      p3m.alpha    = 0.0;
      p3m.alpha_L  = 0.0;
      p3m.r_cut    = 0.0;
      p3m.r_cut_iL = 0.0;
      p3m.mesh[0]  = 0;
      p3m.mesh[1]  = 0;
      p3m.mesh[2]  = 0;
      p3m.cao      = 0;
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
      rf_params.r_cut   = 0.0;
      rf_params.eps   = 0.0;
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

int inter_parse_coulomb(Tcl_Interp * interp, int argc, char ** argv)
{
  double d1;

  Tcl_ResetResult(interp);

  if(argc == 0) {
    printCoulombIAToResult(interp);
    return TCL_OK;
  }
  
  if (! ARG0_IS_D(d1)) {
#ifdef ELP3M
    Tcl_ResetResult(interp);
    if (ARG0_IS_S("elc") && ((coulomb.method == COULOMB_P3M) || (coulomb.method == COULOMB_ELC_P3M)))
      return inter_parse_elc_params(interp, argc - 1, argv + 1);
    if (coulomb.method == COULOMB_P3M)
      return inter_parse_p3m_opt_params(interp, argc, argv);
    else {
      Tcl_AppendResult(interp, "expect: inter coulomb <bjerrum>",
		       (char *) NULL);
      return TCL_ERROR;
    }
#else
    return TCL_ERROR;
#endif
  }

  if (coulomb_set_bjerrum(d1) == TCL_ERROR) {
    Tcl_AppendResult(interp, argv[0], "bjerrum length must be positive",
		     (char *) NULL);
    return TCL_ERROR;
  }
    
  argc -= 1;
  argv += 1;

  if (d1 == 0.0 && argc == 0) {
    mpi_bcast_coulomb_params();
    return TCL_OK;
  }

  if(argc < 1) {
    Tcl_AppendResult(interp, "wrong # args for inter coulomb.",
		     (char *) NULL);
    mpi_bcast_coulomb_params();
    return TCL_ERROR;
  }

  /* check method */

#define REGISTER_COULOMB(name, parser)			\
  if(ARG0_IS_S(name))					\
    return parser(interp, argc-1, argv+1);

#ifdef ELP3M
  REGISTER_COULOMB("p3m", inter_parse_p3m);
#endif

  REGISTER_COULOMB("ewald", inter_parse_ewald);

  REGISTER_COULOMB("dh", inter_parse_dh);    

  REGISTER_COULOMB("rf", inter_parse_rf);

  REGISTER_COULOMB("mmm1d", inter_parse_mmm1d);

  REGISTER_COULOMB("mmm2d", inter_parse_mmm2d);

  REGISTER_COULOMB("maggs", inter_parse_maggs);

  /* fallback */
  coulomb.method  = COULOMB_NONE;
  coulomb.bjerrum = 0.0;

  mpi_bcast_coulomb_params();

  Tcl_AppendResult(interp, "do not know coulomb method \"",argv[0],
		   "\": coulomb switched off", (char *) NULL);
  
  return TCL_ERROR;
}
#endif

/********************************************************************************/
/*                                       printing                               */
/********************************************************************************/

int printBondedIAToResult(Tcl_Interp *interp, int i)
{
  Bonded_ia_parameters *params = &bonded_ia_params[i];
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];

  sprintf(buffer, "%d ", i);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  switch (params->type) {
  case BONDED_IA_FENE:
    Tcl_PrintDouble(interp, params->p.fene.k, buffer);
    Tcl_AppendResult(interp, "FENE ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.fene.r, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
  case BONDED_IA_HARMONIC:
    Tcl_PrintDouble(interp, params->p.harmonic.k, buffer);
    Tcl_AppendResult(interp, "HARMONIC ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.harmonic.r, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
  case BONDED_IA_ANGLE:
    Tcl_PrintDouble(interp, params->p.angle.bend, buffer);
    Tcl_AppendResult(interp, "angle ", buffer," ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.angle.phi0, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
  case BONDED_IA_DIHEDRAL:  
    sprintf(buffer, "%d", params->p.dihedral.mult);
    Tcl_AppendResult(interp, "dihedral ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.dihedral.bend, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.dihedral.phase, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
#ifdef TABULATED
  case BONDED_IA_TABULATED:
    switch (params->p.tab.type) {
    case TAB_BOND_LENGTH:
      Tcl_AppendResult(interp, "tabulated bond \"",params->p.tab.filename,"\"",(char *) NULL);
      return (TCL_OK);
    case TAB_BOND_ANGLE:
      Tcl_AppendResult(interp, "tabulated angle \"",params->p.tab.filename,"\"",(char *) NULL);
      return (TCL_OK);
    case TAB_BOND_DIHEDRAL:
      Tcl_AppendResult(interp, "tabulated dihedral \"",params->p.tab.filename,"\"",(char *) NULL);
      return (TCL_OK);
    }
#endif
#ifdef BOND_CONSTRAINT
  case BONDED_IA_RIGID_BOND:
    Tcl_PrintDouble(interp, sqrt(params->p.rigid_bond.d2), buffer);
    Tcl_AppendResult(interp, "RIGID_BOND ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.rigid_bond.p_tol/2.0, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.rigid_bond.v_tol/time_step, buffer);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return (TCL_OK);
#endif
#ifdef BOND_VIRTUAL
  case BONDED_IA_VIRTUAL_BOND:
    Tcl_AppendResult(interp, "VIRTUAL_BOND ", (char *) NULL);
    return (TCL_OK);
#endif
#ifdef LENNARD_JONES
  case BONDED_IA_SUBT_LJ:
    Tcl_PrintDouble(interp, params->p.subt_lj.k, buffer);
    Tcl_AppendResult(interp, "SUBT_LJ ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, params->p.subt_lj.r, buffer);
    Tcl_AppendResult(interp, buffer,(char *) NULL);
    return (TCL_OK);
#endif
 case BONDED_IA_NONE:
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "unknown bonded interaction number ",buffer,
		     (char *) NULL);
    return (TCL_ERROR);
  }
  /* if none of the above */
  Tcl_ResetResult(interp);
  Tcl_AppendResult(interp, "unknown bonded interaction type",(char *) NULL);
  return (TCL_ERROR);
}

int printNonbondedIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  if (!data) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "interaction does not exist",
		     (char *) NULL);
    return (TCL_ERROR);
  }

  sprintf(buffer, "%d %d ", i, j);
  Tcl_AppendResult(interp, buffer, (char *) NULL);
#ifdef LENNARD_JONES
  if (data->LJ_cut != 0) printljIAToResult(interp,i,j);
#endif
#ifdef LENNARD_JONES_GENERIC
  if (data->LJGEN_cut != 0) printljgenIAToResult(interp,i,j);
#endif
#ifdef SMOOTH_STEP
  if (data->SmSt_cut != 0) printSmStIAToResult(interp,i,j);
#endif
#ifdef BMHTF_NACL
  if (data->BMHTF_cut != 0) printBMHTFIAToResult(interp,i,j);
#endif
#ifdef MORSE
  if (data->MORSE_cut != 0) printmorseIAToResult(interp,i,j);
#endif
#ifdef LJCOS
  if (data->LJCOS_cut != 0) printljcosIAToResult(interp,i,j);
#endif
#ifdef BUCKINGHAM
  if (data->BUCK_cut != 0) printbuckIAToResult(interp,i,j);
#endif
#ifdef SOFT_SPHERE
  if (data->soft_cut != 0) printsoftIAToResult(interp,i,j);
#endif
#ifdef LJCOS2
  if (data->LJCOS2_cut != 0) printljcos2IAToResult(interp,i,j);
#endif
#ifdef ROTATION
  if (data->GB_cut != 0) printgbIAToResult(interp,i,j);
#endif
#ifdef TABULATED
  if (data->TAB_maxval != 0)
    Tcl_AppendResult(interp, "tabulated \"", data->TAB_filename,"\"", (char *) NULL);
#endif
#ifdef COMFORCE
  if (data->COMFORCE_flag != 0) printcomforceIAToResult(interp,i,j);
#endif

#ifdef COMFIXED
  if (data->COMFIXED_flag != 0) printcomfixedIAToResult(interp,i,j);
#endif

#ifdef INTER_DPD
  if ((data->dpd_r_cut != 0)||(data->dpd_tr_cut != 0)) printinterdpdIAToResult(interp,i,j);
#endif

#ifdef INTER_RF
  if (data->rf_r_cut != 0) printinterrfIAToResult(interp,i,j);
#endif
  
  return (TCL_OK);
}

int printCoulombIAToResult(Tcl_Interp *interp) 
{
#ifdef ELECTROSTATICS
  char buffer[TCL_DOUBLE_SPACE + 2*TCL_INTEGER_SPACE];
  if (coulomb.method == COULOMB_NONE) {
    Tcl_AppendResult(interp, "coulomb 0.0", (char *) NULL);
    return (TCL_OK);
  }
  Tcl_PrintDouble(interp, coulomb.bjerrum, buffer);
  Tcl_AppendResult(interp, "{coulomb ", buffer, " ", (char *) NULL);
  switch (coulomb.method) {
#ifdef ELP3M
  case COULOMB_ELC_P3M:
    printP3MToResult(interp);
    printELCToResult(interp);
    break;
  case COULOMB_P3M: printP3MToResult(interp); break;
#endif
  case COULOMB_EWALD: printEWALDToResult(interp); break;
  case COULOMB_DH: printdhToResult(interp); break;
  case COULOMB_RF: printrfToResult(interp); break;
  case COULOMB_MMM1D: printMMM1DToResult(interp); break;
  case COULOMB_MMM2D: printMMM2DToResult(interp); break;
  case COULOMB_MAGGS: printMaggsToResult(interp); break;
  default: break;
  }
  Tcl_AppendResult(interp, "}",(char *) NULL);

#else
  Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)",(char *) NULL);
#endif
  return (TCL_OK);
}

int inter_print_all(Tcl_Interp *interp)
{
  int i, j, start = 1;

  for (i = 0; i < n_bonded_ia; i++) {
    if (bonded_ia_params[i].type != BONDED_IA_NONE) {
      if (start) {
	Tcl_AppendResult(interp, "{", (char *)NULL);
	start = 0;
      }
      else
	Tcl_AppendResult(interp, " {", (char *)NULL);

      printBondedIAToResult(interp, i);
      Tcl_AppendResult(interp, "}", (char *)NULL);
    }
  }
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      if (checkIfParticlesInteract(i, j)) {
	if (start) {
	  Tcl_AppendResult(interp, "{", (char *)NULL);
	  start = 0;
	}
	else
	  Tcl_AppendResult(interp, " {", (char *)NULL);
	printNonbondedIAToResult(interp, i, j);
	Tcl_AppendResult(interp, "}", (char *)NULL);
      }
    }
#ifdef ELECTROSTATICS
  if(coulomb.method != COULOMB_NONE) {
    if (start) 
      start = 0;
    else
      Tcl_AppendResult(interp, " ", (char *)NULL);
    /* here the curled braces will be set inside \ref printCoulombIAToResult
       because electrostatics might be using several lists */
    printCoulombIAToResult(interp);
  }
#endif
  if(lj_force_cap != 0.0) {
    char buffer[TCL_DOUBLE_SPACE];
    
    if (start) {
      Tcl_AppendResult(interp, "{", (char *)NULL);
      start = 0;
    }
    else
      Tcl_AppendResult(interp, " {", (char *)NULL);
    if (lj_force_cap == -1.0)
      Tcl_AppendResult(interp, "ljforcecap individual");
    else {
      Tcl_PrintDouble(interp, lj_force_cap, buffer);
      Tcl_AppendResult(interp, "ljforcecap ", buffer, (char *) NULL);
    }
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }

#ifdef MORSE
if(morse_force_cap != 0.0) {
    char buffer[TCL_DOUBLE_SPACE];

    if (start) {
      Tcl_AppendResult(interp, "{", (char *)NULL);
      start = 0;
    }
    else
      Tcl_AppendResult(interp, " {", (char *)NULL);
    if (morse_force_cap == -1.0)
      Tcl_AppendResult(interp, "morseforcecap individual");
    else {
      Tcl_PrintDouble(interp, morse_force_cap, buffer);
      Tcl_AppendResult(interp, "morseforcecap ", buffer, (char *) NULL);
    }
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
#endif

#ifdef BUCKINGHAM 
  if(buck_force_cap != 0.0) {
    char buffer[TCL_DOUBLE_SPACE];
    
    if (start) {
      Tcl_AppendResult(interp, "{", (char *)NULL);
      start = 0;
    }
    else
      Tcl_AppendResult(interp, " {", (char *)NULL);

    if (lj_force_cap == -1.0)
      Tcl_AppendResult(interp, "buckforcecap individual");
    else {
      Tcl_PrintDouble(interp, buck_force_cap, buffer);
      Tcl_AppendResult(interp, "buckforcecap ", buffer, (char *) NULL);
    }
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
#endif

#ifdef TABULATED
  if(tab_force_cap != 0.0) {
    char buffer[TCL_DOUBLE_SPACE];
    if (start) {
      Tcl_AppendResult(interp, "{", (char *)NULL);
      start = 0;
    }
    else
      Tcl_AppendResult(interp, " {", (char *)NULL);
    if (tab_force_cap == -1.0)
      Tcl_AppendResult(interp, "tabforcecap individual");
    else {
      Tcl_PrintDouble(interp, tab_force_cap, buffer);
      Tcl_AppendResult(interp, "tabforcecap ", buffer, (char *) NULL);
    }
    Tcl_AppendResult(interp, "}", (char *)NULL);
  }
#endif

  return (TCL_OK);
}

int inter_print_bonded(Tcl_Interp *interp, int i)
{
  char buffer[TCL_INTEGER_SPACE];

  Tcl_ResetResult(interp);

  if(i < 0) {
    Tcl_AppendResult(interp, "interaction type must be nonnegative",
		     (char *) NULL);
    return (TCL_ERROR);
  }
  
  /* print specific interaction information */
  if(i<n_bonded_ia) {
    printBondedIAToResult(interp, i);
    return TCL_OK;
  }

  sprintf(buffer, "%d", i);
  Tcl_AppendResult(interp, "unknown bonded interaction number ", buffer,
		   (char *) NULL);
  return TCL_ERROR;
}

int inter_print_non_bonded(Tcl_Interp * interp,
			   int part_type_a, int part_type_b)
{
  IA_parameters *data, *data_sym;

  Tcl_ResetResult(interp);

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    Tcl_AppendResult(interp, "particle types must be nonnegative",
		     (char *) NULL);
    return TCL_ERROR;
  }

  return printNonbondedIAToResult(interp, part_type_a, part_type_b);
}

int inter_parse_non_bonded(Tcl_Interp * interp,
			   int part_type_a, int part_type_b,
			   int argc, char ** argv)
{
  int change;
  
  Tcl_ResetResult(interp);

  if (argc <= 0) {
    Tcl_AppendResult(interp, "wrong # args:  should be \"",
		     "inter <type 1> <type 2> ?interaction? ?values?\"",
		     (char *) NULL);
    return TCL_ERROR;
  }

  /* get interaction parameters */

  while (argc > 0) {
    /* The various parsers return the number of parsed parameters.
       If an error occured, 0 should be returned, since none of the parameters were
       understood */

    /* that's just for the else below... */
    if (0);

#define REGISTER_NONBONDED(name, parser)				\
    else if (ARG0_IS_S(name))						\
      change = parser(interp, part_type_a, part_type_b, argc, argv)

#ifdef LENNARD_JONES
    REGISTER_NONBONDED("lennard-jones", lj_parser);
#endif

#ifdef LENNARD_JONES_GENERIC
    REGISTER_NONBONDED("lj-gen", ljgen_parser);
#endif

#ifdef SMOOTH_STEP
    REGISTER_NONBONDED("smooth-step", SmSt_parser);
#endif

#ifdef BMHTF_NACL
    REGISTER_NONBONDED("bmhtf-nacl", BMHTF_parser);
#endif

#ifdef MORSE
    REGISTER_NONBONDED("morse", morse_parser);
#endif

#ifdef LJCOS
    REGISTER_NONBONDED("lj-cos", ljcos_parser);
#endif

#ifdef BUCKINGHAM
    REGISTER_NONBONDED("buckingham", buckingham_parser);
#endif

#ifdef SOFT_SPHERE
    REGISTER_NONBONDED("soft-sphere", soft_parser);
#endif

#ifdef COMFORCE
    REGISTER_NONBONDED("comforce", comforce_parser);
#endif

#ifdef LJCOS2
    REGISTER_NONBONDED("lj-cos2", ljcos2_parser);
#endif

#ifdef COMFIXED
    REGISTER_NONBONDED("comfixed", comfixed_parser);
#endif

#ifdef ROTATION
    REGISTER_NONBONDED("gay-berne", gb_parser);
#endif

#ifdef TABULATED
    REGISTER_NONBONDED("tabulated", tab_parser);
#endif
#ifdef INTER_DPD
    REGISTER_NONBONDED("inter_dpd", interdpd_parser);
#endif
#ifdef INTER_RF
    REGISTER_NONBONDED("inter_rf", interrf_parser);
#endif
    else {
      Tcl_AppendResult(interp, "excessive parameter/unknown interaction type \"", argv[0],
		       "\" in parsing non bonded interaction",
		       (char *) NULL);
      return TCL_ERROR;
    }

    if (change <= 0)
      return TCL_ERROR;

    argc -= change;
    argv += change;
  }
  return TCL_OK;
}

int inter_print_partner_num(Tcl_Interp *interp, int bond_type)
{
  Bonded_ia_parameters * params;
  char buffer[TCL_INTEGER_SPACE];

  if(bond_type < 0) {
    Tcl_AppendResult(interp, "interaction type must be nonnegative",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if(bond_type < n_bonded_ia) {
    params = &bonded_ia_params[bond_type];
    sprintf(buffer, "%d", params->num);
    Tcl_AppendResult(interp, buffer, (char *) NULL);
    return TCL_OK;
  }
 
  sprintf(buffer, "%d", bond_type);
  Tcl_AppendResult(interp, "unknown bonded interaction number ", buffer,
		   (char *) NULL);
  return TCL_ERROR;
}

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

int inter_parse_virtual_bonds(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
	CHECK_VALUE(virtual_set_params(bond_type), "bond type must be nonnegative");
}
#endif

int inter_parse_bonded(Tcl_Interp *interp,
		       int bond_type,
		       int argc, char ** argv)
{
  if (ARG0_IS_S("num")) {
    if (argc == 1)
      return inter_print_partner_num(interp, bond_type);
    else {
	Tcl_AppendResult(interp, "to manny parameters",
			 (char *) NULL);
	return TCL_ERROR;
    }
  }

#define REGISTER_BONDED(name,parser)			\
  if (ARG0_IS_S(name)) return parser(interp, bond_type, argc, argv);
  
  REGISTER_BONDED("fene", inter_parse_fene);
  REGISTER_BONDED("harmonic", inter_parse_harmonic);
#ifdef LENNARD_JONES  
  REGISTER_BONDED("subt_lj", inter_parse_subt_lj);
#endif
#ifdef BOND_ANGLE
  REGISTER_BONDED("angle", inter_parse_angle);
#endif
  REGISTER_BONDED("dihedral", inter_parse_dihedral);
#ifdef TABULATED
  REGISTER_BONDED("tabulated", inter_parse_bonded_tabulated);
#endif
#ifdef BOND_CONSTRAINT
  REGISTER_BONDED("rigid_bond", inter_parse_rigid_bonds);
#endif
#ifdef BOND_VIRTUAL
  REGISTER_BONDED("virtual_bond", inter_parse_virtual_bonds);
#endif
  Tcl_AppendResult(interp, "unknown bonded interaction type \"", argv[0],
		   "\"", (char *) NULL);
  return TCL_ERROR;
}

int inter_parse_rest(Tcl_Interp * interp, int argc, char ** argv)
{
#if defined(LENNARD_JONES) || defined(LENNARD_JONES_GENERIC)
  if(ARG0_IS_S("ljforcecap"))
    return inter_parse_ljforcecap(interp, argc-1, argv+1);
#endif

#ifdef MORSE
  if(ARG0_IS_S("morseforcecap"))
    return inter_parse_morseforcecap(interp, argc-1, argv+1);
#endif

#ifdef BUCKINGHAM
  if(ARG0_IS_S("buckforcecap"))
    return inter_parse_buckforcecap(interp, argc-1, argv+1);
#endif

#ifdef TABULATED
  if(ARG0_IS_S("tabforcecap"))
    return inter_parse_tabforcecap(interp, argc-1, argv+1);
#endif

  if(ARG0_IS_S("coulomb")) {
#ifdef ELECTROSTATICS
    return inter_parse_coulomb(interp, argc-1, argv+1);
#else
    Tcl_AppendResult(interp, "ELECTROSTATICS not compiled (see config.h)", (char *) NULL);
#endif
  }
  Tcl_AppendResult(interp, "unknown interaction type \"", argv[0],
		   "\"", (char *) NULL);

  return TCL_ERROR;
}

int inter(ClientData _data, Tcl_Interp *interp,
	  int argc, char **argv)
{
  int i, j, err_code = TCL_OK, is_i1, is_i2;

  Tcl_ResetResult(interp);

  /* first we handle the special cases

     1. no parameters
     2. one parameter
     3. two parameters

     then the bonded interactions
     then the non bonded interactions
     then the rest
   */

  if (argc == 1) {
    /* no argument -> print all interaction informations. */
    err_code = inter_print_all(interp);
  }
  else if (argc == 2) {
    /* There is only 1 parameter, bonded ia printing or force caps */

    if (ARG1_IS_I(i))
      err_code = inter_print_bonded(interp, i);
    else {
      Tcl_ResetResult(interp);
      err_code = inter_parse_rest(interp, argc-1, argv+1);
    }
  }
  else if (argc == 3) {
    /* There are only 2 parameters, non_bonded printing */
    
    is_i1 = ARG_IS_I(1, i);
    is_i2 = ARG_IS_I(2, j);

    Tcl_ResetResult(interp);

    if (is_i1 && is_i2)
      err_code = inter_print_non_bonded(interp, i, j);
    else if (is_i1)
      err_code = inter_parse_bonded(interp, i, argc-2, argv+2);
    else
      err_code = inter_parse_rest(interp, argc-1, argv+1);
  }
  else {
    /****************************************************
     * Here we have more than 2 parameters
     ****************************************************/

    is_i1 = ARG_IS_I(1, i);
    is_i2 = ARG_IS_I(2, j);

    Tcl_ResetResult(interp);
 
    // non bonded interactions
    if (is_i1 && is_i2)
      err_code = inter_parse_non_bonded(interp, i, j, argc-3, argv+3);
    else if (is_i1)
      // bonded interactions
      err_code = inter_parse_bonded(interp, i, argc-2, argv+2);
    else
      // named interactions
      err_code = inter_parse_rest(interp, argc-1, argv+1);
  }
  /* check for background errors which have not been handled so far */
  return mpi_gather_runtime_errors(interp, err_code);
}
