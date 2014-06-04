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
/** \file interaction_data.cpp
    Implementation of interaction_data.hpp
 */
#include <cstring>
#include <cstdlib>
#include "utils.hpp"
#include "rattle.hpp"
#include "interaction_data.hpp"
#include "errorhandling.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "pressure.hpp"
#include "p3m.hpp"
#include "debye_hueckel.hpp"
#include "reaction_field.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "maggs.hpp"
#include "elc.hpp"
#include "lj.hpp"
#include "ljgen.hpp"
#include "ljangle.hpp"
#include "steppot.hpp"
#include "hertzian.hpp"
#include "gaussian.hpp"
#include "buckingham.hpp"
#include "soft_sphere.hpp"
#include "hat.hpp"
#include "tab.hpp"
#include "overlap.hpp"
#include "ljcos.hpp"
#include "ljcos2.hpp"
#include "gb.hpp"
#include "cells.hpp"
#include "comforce.hpp"
#include "comfixed.hpp"
#include "morse.hpp"
#include "dpd.hpp"
#include "tunable_slip.hpp"
#include "magnetic_non_p3m_methods.hpp"
#include "mdlc_correction.hpp"
#include "initialize.hpp"

/****************************************
 * variables
 *****************************************/
int n_particle_types = 0;
int n_interaction_types = 0;
IA_parameters *ia_params = NULL;

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

double min_global_cut = 0.0;

double max_cut;
double max_cut_nonbonded;
double max_cut_bonded;
/** maximal cutoff of type-independent short range ia, mainly
    electrostatics and DPD*/
double max_cut_global;

/** Array containing all tabulated forces*/
DoubleList tabulated_forces;
/** Corresponding array containing all tabulated energies*/
DoubleList tabulated_energies;

/*****************************************
 * function prototypes
 *****************************************/

/** calculates and returns the maximal global nonbonded cutoff that is
    required.  Currently, this are just the cutoffs from the
    electrostatics method and some dpd cutoffs. */
static void recalc_global_maximal_nonbonded_cutoff();
/** calculate the maximal cutoff of bonded interactions, required to
    determine the cell size for communication. */
static void recalc_maximal_cutoff_bonded();

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
 
  params->particlesInteract = 0;
  params->max_cut = max_cut_global;

#ifdef LENNARD_JONES
  params->LJ_eps =
    params->LJ_sig =
    params->LJ_shift =
    params->LJ_offset =
    params->LJ_capradius =
    params->LJ_min = 0.0;
  params->LJ_cut = INACTIVE_CUTOFF;
#endif

#ifdef LENNARD_JONES_GENERIC
  params->LJGEN_eps =
    params->LJGEN_sig =
    params->LJGEN_shift =
    params->LJGEN_offset =
    params->LJGEN_capradius =
    params->LJGEN_a1 =
    params->LJGEN_a2 = 
    params->LJGEN_b1 =
    params->LJGEN_b2 = 0.0;
#ifdef LJGEN_SOFTCORE
    params->LJGEN_lambda  = 1.0;
    params->LJGEN_softrad = 0.0;
#endif
  params->LJGEN_cut = INACTIVE_CUTOFF;
#endif

#ifdef LJ_ANGLE
  params->LJANGLE_eps =
    params->LJANGLE_sig =
    params->LJANGLE_bonded1type=
    params->LJANGLE_bonded1pos = 
    params->LJANGLE_bonded1neg = 
    params->LJANGLE_bonded2pos = 
    params->LJANGLE_bonded2neg = 
    params->LJANGLE_capradius =
    params->LJANGLE_z0 =
    params->LJANGLE_kappa =
    params->LJANGLE_epsprime = 0.0;
  params->LJANGLE_dz = -1.0;
  params->LJANGLE_cut = INACTIVE_CUTOFF;
#endif

#ifdef SMOOTH_STEP
  params->SmSt_eps =
    params->SmSt_sig =
    params->SmSt_d =
    params->SmSt_n =
    params->SmSt_k0 = 0.0;
  params->SmSt_cut = INACTIVE_CUTOFF;
#endif

#ifdef HERTZIAN
  params->Hertzian_eps = 0.0;
  params->Hertzian_sig = INACTIVE_CUTOFF;
#endif

#ifdef GAUSSIAN
  params->Gaussian_eps = 0.0;
  params->Gaussian_sig = 1.0;
  params->Gaussian_cut = INACTIVE_CUTOFF;
#endif

#ifdef BMHTF_NACL
  params->BMHTF_A =
    params->BMHTF_B =
    params->BMHTF_C =
    params->BMHTF_D =
    params->BMHTF_sig =
    params->BMHTF_computed_shift = 0.0;
  params->BMHTF_cut = INACTIVE_CUTOFF;
#endif

#ifdef MORSE
  params->MORSE_eps = 
    params->MORSE_alpha =
    params->MORSE_rmin =
    params->MORSE_rest = 
    params->MORSE_capradius = 0;
  params->MORSE_cut = INACTIVE_CUTOFF;
#endif

#ifdef BUCKINGHAM
  params->BUCK_A =
    params->BUCK_B =
    params->BUCK_C =
    params->BUCK_D =
    params->BUCK_discont =
    params->BUCK_shift =
    params->BUCK_capradius =
    params->BUCK_F1 =
    params->BUCK_F2 = 0.0;
  params->BUCK_cut = INACTIVE_CUTOFF;
#endif

#ifdef SOFT_SPHERE
  params->soft_a =
    params->soft_n =
    params->soft_offset = 0.0;
  params->soft_cut = INACTIVE_CUTOFF;
#endif

#ifdef HAT
  params->HAT_Fmax =
    params->HAT_r = 0.0;
#endif

#ifdef LJCOS
  params->LJCOS_eps =
    params->LJCOS_sig =
    params->LJCOS_offset =
    params->LJCOS_alfa =
    params->LJCOS_beta =
    params->LJCOS_rmin = 0.0;
  params->LJCOS_cut = INACTIVE_CUTOFF;
#endif

#ifdef LJCOS2
  params->LJCOS2_eps =
    params->LJCOS2_sig =
    params->LJCOS2_offset =
    params->LJCOS2_w =
    params->LJCOS2_rchange = 
    params->LJCOS2_capradius = 0.0;
  params->LJCOS2_cut = INACTIVE_CUTOFF;
#endif

#ifdef GAY_BERNE
  params->GB_eps =
    params->GB_sig =
    params->GB_k1 =
    params->GB_k2 =
    params->GB_mu =
    params->GB_nu =
    params->GB_chi1 =
    params->GB_chi2 = 0.0;
  params->GB_cut = INACTIVE_CUTOFF;
#endif

#ifdef TABULATED
  params->TAB_npoints =
    params->TAB_startindex = 0;
  params->TAB_minval =
    params->TAB_stepsize = 0.0;
  strcpy(params->TAB_filename,"");
  params->TAB_maxval = INACTIVE_CUTOFF;
#endif

#ifdef INTER_DPD
  params->dpd_gamma = 0.0;
  params->dpd_wf = 0;
  params->dpd_pref1 = 0.0;
  params->dpd_pref2 = 0.0;
  params->dpd_tgamma = 0.0;
  params->dpd_tr_cut = 0.0;
  params->dpd_wf = 0;
  params->dpd_pref3 = 0;
  params->dpd_pref4 = 0;
  params->dpd_r_cut = INACTIVE_CUTOFF;
#endif

#ifdef TUNABLE_SLIP
  params->TUNABLE_SLIP_temp  = 0.0;
  params->TUNABLE_SLIP_gamma = 0.0;
  params->TUNABLE_SLIP_time  = 0.0;
  params->TUNABLE_SLIP_vx  = 0.0;
  params->TUNABLE_SLIP_vy  = 0.0;
  params->TUNABLE_SLIP_vz  = 0.0;
  params->TUNABLE_SLIP_r_cut = INACTIVE_CUTOFF;
#endif

  /* things that are not strictly speaking short-ranged interactions,
     and do not have a cutoff */
#ifdef COMFORCE
  params->COMFORCE_flag = 0;
  params->COMFORCE_dir = 0;
  params->COMFORCE_force = 0.;
  params->COMFORCE_fratio = 0.;
#endif

#ifdef COMFIXED
  params->COMFIXED_flag = 0;
#endif

#ifdef INTER_RF
  params->rf_on = 0;
#endif

#ifdef MOL_CUT
  params->mol_cut_type = 0;
  params->mol_cut_cutoff = 0.0;
#endif

#ifdef CATALYTIC_REACTIONS
  params->REACTION_range = 0.0;
#endif

}

/** Copy interaction parameters. */
void copy_ia_params(IA_parameters *dst, IA_parameters *src) {
  memcpy(dst, src, sizeof(IA_parameters));
}

IA_parameters *get_ia_param_safe(int i, int j) {
  make_particle_type_exist(imax(i, j));
  return get_ia_param(i, j);
}

static void recalc_maximal_cutoff_bonded()
{
  int i;
  double max_cut_tmp;

  max_cut_bonded = 0.0;

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
}

static void recalc_global_maximal_nonbonded_cutoff()
{
  /* user defined minimal global cut. This makes sure that data of
   pairs of particles with a distance smaller than this are always
   available on the same node (through ghosts). Required for example
   for the relative virtual sites algorithm. */
  max_cut_global = min_global_cut;

#ifdef ELECTROSTATICS
  /* Cutoff for the real space electrostatics.
     Note that the box length may have changed,
     but the method not yet reinitialized.
   */
  switch (coulomb.method) {
#ifdef P3M 
  case COULOMB_ELC_P3M:
    if (max_cut_global < elc_params.space_layer)
      max_cut_global = elc_params.space_layer;
    // fall through
  case COULOMB_P3M_GPU:
  case COULOMB_P3M: {
    /* do not use precalculated r_cut here, might not be set yet */
    double r_cut = p3m.params.r_cut_iL* box_l[0];
    if (max_cut_global < r_cut)
      max_cut_global = r_cut;
    break;
  }
#endif
  case COULOMB_DH:
    if (max_cut_global < dh_params.r_cut)
      max_cut_global = dh_params.r_cut;
    break;
  case COULOMB_RF:
  case COULOMB_INTER_RF:
    if (max_cut_global < rf_params.r_cut)
      max_cut_global = rf_params.r_cut;
    break;
  }
#endif /*ifdef ELECTROSTATICS */
  
#ifdef DIPOLES
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case DIPOLAR_MDLC_P3M:
    // fall through
  case DIPOLAR_P3M: {
    /* do not use precalculated r_cut here, might not be set yet */
    double r_cut = dp3m.params.r_cut_iL* box_l[0];
    if (max_cut_global < r_cut)
      max_cut_global = r_cut;
    break;
  }
#endif /*ifdef DP3M */
  }       
#endif

#ifdef DPD
  if (dpd_r_cut != 0) {
    if(max_cut_global < dpd_r_cut)
      max_cut_global = dpd_r_cut;
  }
#endif
  
#ifdef TRANS_DPD
  if (dpd_tr_cut != 0) {
    if(max_cut_global < dpd_tr_cut)
      max_cut_global = dpd_tr_cut;
  }
#endif

}

static void recalc_maximal_cutoff_nonbonded()
{
  int i, j;

  CELL_TRACE(fprintf(stderr, "%d: recalc_maximal_cutoff_nonbonded\n", this_node));

  recalc_global_maximal_nonbonded_cutoff();

  CELL_TRACE(fprintf(stderr, "%d: recalc_maximal_cutoff_nonbonded: max_cut_global = %f\n", this_node, max_cut_global));

  max_cut_nonbonded = max_cut_global;
  
  for (i = 0; i < n_particle_types; i++)
    for (j = i; j < n_particle_types; j++) {
      double max_cut_current = 0;

      IA_parameters *data = get_ia_param(i, j);

#ifdef LENNARD_JONES
      if(max_cut_current < (data->LJ_cut+data->LJ_offset))
	max_cut_current = (data->LJ_cut+data->LJ_offset);
#endif

#ifdef INTER_DPD
      {
	double max_cut_tmp = (data->dpd_r_cut > data->dpd_tr_cut) ?
	  data->dpd_r_cut : data->dpd_tr_cut;
	if (max_cut_current <  max_cut_tmp)
	  max_cut_current = max_cut_tmp;
      }
#endif

#ifdef LENNARD_JONES_GENERIC
      if (max_cut_current < (data->LJGEN_cut+data->LJGEN_offset))
	max_cut_current = (data->LJGEN_cut+data->LJGEN_offset);
#endif

#ifdef LJ_ANGLE
      if (max_cut_current < (data->LJANGLE_cut))
	max_cut_current = (data->LJANGLE_cut);
#endif

#ifdef SMOOTH_STEP
      if (max_cut_current < data->SmSt_cut)
	max_cut_current = data->SmSt_cut;
#endif

#ifdef HERTZIAN
      if (max_cut_current < data->Hertzian_sig)
	max_cut_current = data->Hertzian_sig;
#endif

#ifdef GAUSSIAN
      if (max_cut_current < data->Gaussian_cut)
	max_cut_current = data->Gaussian_cut;
#endif

#ifdef BMHTF_NACL
      if (max_cut_current < data->BMHTF_cut)
	max_cut_current = data->BMHTF_cut;
#endif

#ifdef MORSE
      if (max_cut_current < data->MORSE_cut)
	max_cut_current = data->MORSE_cut;
#endif

#ifdef BUCKINGHAM
      if (max_cut_current < data->BUCK_cut)
	max_cut_current = data->BUCK_cut;
#endif

#ifdef SOFT_SPHERE
      if (max_cut_current < data->soft_cut)
	max_cut_current = data->soft_cut;
#endif

#ifdef HAT
      if (max_cut_current < data->HAT_r)
	max_cut_current = data->HAT_r;
#endif

#ifdef LJCOS
      {
	double max_cut_tmp = data->LJCOS_cut + data->LJCOS_offset;
	if (max_cut_current < max_cut_tmp)
	  max_cut_current = max_cut_tmp;
      }
#endif

#ifdef LJCOS2
      {
	double max_cut_tmp = data->LJCOS2_cut + data->LJCOS2_offset;
	if (max_cut_current < max_cut_tmp)
	  max_cut_current = max_cut_tmp;
      }
#endif

#ifdef GAY_BERNE
      if (max_cut_current < data->GB_cut)
	max_cut_current = data->GB_cut;
#endif

#ifdef TABULATED
      if (max_cut_current < data->TAB_maxval)
	max_cut_current = data->TAB_maxval;
#endif
	 
#ifdef TUNABLE_SLIP
      if (max_cut_current < data->TUNABLE_SLIP_r_cut)
	max_cut_current = data->TUNABLE_SLIP_r_cut;
#endif

#ifdef CATALYTIC_REACTIONS
      if (max_cut_current < data->REACTION_range)
	max_cut_current = data->REACTION_range;
#endif

#ifdef MOL_CUT
      if (data->mol_cut_type != 0) {
	if (max_cut_current < data->mol_cut_cutoff)
	  max_cut_current = data->mol_cut_cutoff;
	max_cut_current += 2.0* max_cut_bonded;
      }
#endif

      IA_parameters *data_sym = get_ia_param(j, i);

      /* no interaction ever touched it, at least no real
	 short-ranged one (that writes to the nonbonded energy) */
      data_sym->particlesInteract =
	data->particlesInteract = (max_cut_current > 0.0);
      
      /* take into account any electrostatics */
      if (max_cut_global > max_cut_current)
	max_cut_current = max_cut_global;

      data_sym->max_cut =
	data->max_cut = max_cut_current;

      if (max_cut_current > max_cut_nonbonded)
	max_cut_nonbonded = max_cut_current;

      CELL_TRACE(fprintf(stderr, "%d: pair %d,%d max_cut total %f\n",
			 this_node, i, j, data->max_cut));
    }
}

void recalc_maximal_cutoff()
{
  recalc_maximal_cutoff_bonded();
  recalc_maximal_cutoff_nonbonded();

  /* make max_cut the maximal cutoff of both bonded and non-bonded
     interactions */
  if (max_cut_nonbonded > max_cut_bonded)
    max_cut = max_cut_nonbonded;
  else
    max_cut = max_cut_bonded;
}

const char *get_name_of_bonded_ia(int i) {
  switch (i) {
  case BONDED_IA_FENE:
    return "FENE";
  case BONDED_IA_ANGLE_OLD:
    return "angle";
  case BONDED_IA_ANGLE_HARMONIC:
    return "angle_harmonic";
  case BONDED_IA_ANGLE_COSINE:
    return "angle_cosine";
  case BONDED_IA_ANGLE_COSSQUARE:
    return "angle_cossquare";
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
  case BONDED_IA_STRETCHING_FORCE:
    return "STRETCHING_FORCE";
  case BONDED_IA_AREA_FORCE_LOCAL:
    return "AREA_FORCE_LOCAL";
  case BONDED_IA_AREA_FORCE_GLOBAL:
    return "AREA_FORCE_GLOBAL";
  case BONDED_IA_BENDING_FORCE:
    return "BENDING_FORCE";
  case BONDED_IA_VOLUME_FORCE:
    return "VOLUME_FORCE";
  case BONDED_IA_STRETCHLIN_FORCE:
    return "STRETCHLIN_FORCE";
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

int interactions_sanity_checks()
{
  /* set to zero if initialization was not successful. */
  int state = 1;

#ifdef ELECTROSTATICS
  switch (coulomb.method) {
  case COULOMB_MMM1D: if (MMM1D_sanity_checks()) state = 0; break;
  case COULOMB_MMM2D: if (MMM2D_sanity_checks()) state = 0; break;
#ifdef P3M
  case COULOMB_ELC_P3M: if (ELC_sanity_checks()) state = 0; // fall through
  case COULOMB_P3M_GPU:
  case COULOMB_P3M: if (p3m_sanity_checks()) state = 0; break;
#endif
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
    return ES_ERROR;
  
  coulomb.bjerrum = bjerrum;

  if (coulomb.bjerrum == 0.0) {
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_ELC_P3M:
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
      p3m_set_bjerrum();
      break;
#endif
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
    }
 
    mpi_bcast_coulomb_params();
    coulomb.method = COULOMB_NONE;
    mpi_bcast_coulomb_params();

  }

  return ES_OK;
}



/* =========================================================
   ========================================================= */
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES

int dipolar_set_Dbjerrum(double bjerrum)
{
  if (bjerrum < 0.0)
    return ES_ERROR;
  
  coulomb.Dbjerrum = bjerrum;

  if (coulomb.Dbjerrum == 0.0) {
    switch (coulomb.Dmethod) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:
      // fall through
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

  return ES_OK;
}

#endif   /* ifdef  DIPOLES */

void recalc_coulomb_prefactor()
{
#ifdef ELECTROSTATICS
  if(temperature > 0.0)
    coulomb.prefactor = coulomb.bjerrum * temperature; 
  else
    coulomb.prefactor = coulomb.bjerrum;
#endif

#ifdef DIPOLES
  if(temperature > 0.0)
    coulomb.Dprefactor = coulomb.Dbjerrum * temperature; 
  else
    coulomb.Dprefactor = coulomb.Dbjerrum;
#endif
}

#ifdef BOND_VIRTUAL
int virtual_set_params(int bond_type)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].type = BONDED_IA_VIRTUAL_BOND;
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

#endif
