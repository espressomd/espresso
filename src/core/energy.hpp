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
/** \file energy.hpp
    Implementation of the energy calculation.
*/
#ifndef _ENERGY_H
#define _ENERGY_H
#include "utils.hpp"
#include "integrate.hpp"
#include "statistics.hpp"
#include "thermostat.hpp"

/* include the energy files */
#include "p3m.hpp"
#include "p3m-dipolar.hpp"
#include "lj.hpp"
#include "ljgen.hpp"
#include "steppot.hpp"
#include "hertzian.hpp"
#include "gaussian.hpp"
#include "bmhtf-nacl.hpp"
#include "buckingham.hpp"
#include "soft_sphere.hpp"
#include "hat.hpp"
#include "ljcos.hpp"
#include "ljcos2.hpp"
#include "ljangle.hpp"
#include "tab.hpp"
#include "overlap.hpp"
#include "gb.hpp"
#include "fene.hpp"
#include "object-in-fluid/stretching_force.hpp"
#include "object-in-fluid/stretchlin_force.hpp"
#include "object-in-fluid/area_force_local.hpp"
#include "object-in-fluid/area_force_global.hpp"
#include "object-in-fluid/bending_force.hpp"
#include "object-in-fluid/volume_force.hpp"
#include "harmonic.hpp"
#include "subt_lj.hpp"
#include "angle.hpp"
#include "angle_harmonic.hpp"
#include "angle_cosine.hpp"
#include "angle_cossquare.hpp"
#include "angledist.hpp"
#include "dihedral.hpp"
#include "debye_hueckel.hpp"
#include "endangledist.hpp"
#include "reaction_field.hpp"
#include "mmm1d.hpp"
#include "mmm2d.hpp"
#include "maggs.hpp"
#include "morse.hpp"
#include "elc.hpp"
#include "mdlc_correction.hpp"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat energy, total_energy;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** allocate energy arrays and initialize with zero */
void init_energies(Observable_stat *stat);

/** on the master node: calc energies only if necessary */
void master_energy_calc();

/** parallel energy calculation.
    @param result non-zero only on master node; will contain the cumulative over all nodes. */
void energy_calc(double *result);

/** Calculate non bonded energies between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param ia_params the interaction parameters between the two particles
    @param d         vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2.
    @return the short ranged interaction energy between the two particles
*/
inline double calc_non_bonded_pair_energy(Particle *p1, Particle *p2,
					    IA_parameters *ia_params,
					    double d[3], double dist, double dist2)
{
  double ret = 0;

#ifdef NO_INTRA_NB
  if (p1->p.mol_id==p2->p.mol_id) return 0;
#endif

#ifdef MOL_CUT
  //You may want to put a correction factor for smoothing function else then theta
  if (checkIfParticlesInteractViaMolCut(p1,p2,ia_params)==0) return 0;
#endif

#ifdef LENNARD_JONES
  /* lennard jones */
  ret += lj_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef LENNARD_JONES_GENERIC
  /* Generic lennard jones */
  ret += ljgen_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef LJ_ANGLE
  /* Directional LJ */
  ret += ljangle_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef SMOOTH_STEP
  /* smooth step */
  ret += SmSt_pair_energy(p1,p2,ia_params,d,dist,dist2);
#endif

#ifdef HERTZIAN
  /* Hertzian potential */
  ret += hertzian_pair_energy(p1,p2,ia_params,d,dist,dist2);
#endif

#ifdef GAUSSIAN
  /* Gaussian potential */
  ret += gaussian_pair_energy(p1,p2,ia_params,d,dist,dist2);
#endif

#ifdef BMHTF_NACL
  /* BMHTF NaCl */
  ret += BMHTF_pair_energy(p1,p2,ia_params,d,dist,dist2);
#endif

#ifdef MORSE
  /* morse */
  ret +=morse_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef BUCKINGHAM
  /* lennard jones */
  ret  += buck_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef SOFT_SPHERE
  /* soft-sphere */
  ret  += soft_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef HAT
  /* hat */
  ret  += hat_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef LJCOS2
  /* lennard jones */
  ret += ljcos2_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef TABULATED
  /* tabulated */
  ret += tabulated_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef LJCOS
  /* lennard jones cosine */
  ret += ljcos_pair_energy(p1,p2,ia_params,d,dist);
#endif
  
#ifdef GAY_BERNE
  /* Gay-Berne */
  ret += gb_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef INTER_RF
  ret += interrf_pair_energy(p1,p2,ia_params,dist);
#endif

  return ret;
}

/** Add non bonded energies and short range coulomb between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
inline void add_non_bonded_pair_energy(Particle *p1, Particle *p2, double d[3],
					 double dist, double dist2)
{
  IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);

#if defined(ELECTROSTATICS) || defined(DIPOLES)
  double ret = 0;
#endif

  *obsstat_nonbonded(&energy, p1->p.type, p2->p.type) +=
    calc_non_bonded_pair_energy(p1, p2, ia_params, d, dist, dist2);

#ifdef ELECTROSTATICS
  if (coulomb.method != COULOMB_NONE) {
    /* real space coulomb */
    switch (coulomb.method) {
#ifdef P3M
    case COULOMB_P3M_GPU:
    case COULOMB_P3M:
      ret = p3m_pair_energy(p1->p.q*p2->p.q,d,dist2,dist);
      break;
    case COULOMB_ELC_P3M:
      ret = p3m_pair_energy(p1->p.q*p2->p.q,d,dist2,dist);
      if (elc_params.dielectric_contrast_on)
      ret += 0.5*ELC_P3M_dielectric_layers_energy_contribution(p1,p2);
    break;
#endif
    case COULOMB_DH:
      ret = dh_coulomb_pair_energy(p1,p2,dist);
      break;
    case COULOMB_RF:
      ret = rf_coulomb_pair_energy(p1,p2,dist);
      break;
    case COULOMB_INTER_RF:
      //this is done above as interaction
      ret = 0;
      break;
    case COULOMB_MMM1D:
      ret = mmm1d_coulomb_pair_energy(p1,p2,d,dist2,dist);
      break;
    case COULOMB_MMM2D:
      ret = mmm2d_coulomb_pair_energy(p1->p.q*p2->p.q,d,dist2,dist);
      break;
    default :
      ret = 0.;
    }
    energy.coulomb[0] += ret;
  }
#endif

#ifdef DIPOLES
  if (coulomb.Dmethod != DIPOLAR_NONE) {
    ret=0;
    switch (coulomb.Dmethod) {
#ifdef DP3M
    case DIPOLAR_MDLC_P3M:  
      //fall trough
    case DIPOLAR_P3M:
      ret = dp3m_pair_energy(p1,p2,d,dist2,dist); 
      break;
#endif 
    }
    energy.dipolar[0] += ret;
  }
#endif

}

/** Calculate bonded energies for one particle.
    @param p1 particle for which to calculate energies
*/
inline void add_bonded_energy(Particle *p1)
{
  char *errtxt;
  Particle *p2, *p3 = NULL, *p4 = NULL;
  Bonded_ia_parameters *iaparams;
  int i, type_num, type, n_partners, bond_broken;
  double ret=0, dx[3] = {0, 0, 0};

  i = 0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];
    type = iaparams->type;
    n_partners = iaparams->num;
    
    /* fetch particle 2, which is always needed */
    p2 = local_particles[p1->bl.e[i++]];
    if (!p2) {
      errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{069 bond broken between particles %d and %d (particles not stored on the same node)} ",
	      p1->p.identity, p1->bl.e[i-1]);
      return;
    }

    /* fetch particle 3 eventually */
    if (n_partners >= 2) {
      p3 = local_particles[p1->bl.e[i++]];
      if (!p3) {
	errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{070 bond broken between particles %d, %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i-2], p1->bl.e[i-1]);
	return;
      }
    }

    /* fetch particle 4 eventually */
    if (n_partners >= 3) {
      p4 = local_particles[p1->bl.e[i++]];
      if (!p4) {
	errtxt = runtime_error(128 + 4*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{071 bond broken between particles %d, %d, %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i-3], p1->bl.e[i-2], p1->bl.e[i-1]);
	return;
      }
    }

    /* similar to the force, we prepare the center-center vector */
    if (n_partners == 1)
      get_mi_vector(dx, p1->r.p, p2->r.p);

    switch(type) {
    case BONDED_IA_FENE:
      bond_broken = fene_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
    case BONDED_IA_HARMONIC:
      bond_broken = harmonic_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#ifdef LENNARD_JONES
    case BONDED_IA_SUBT_LJ:
      bond_broken = subt_lj_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#endif
#ifdef BOND_ANGLE_OLD
    /* the first case is not needed and should not be called */ 
    case BONDED_IA_ANGLE_OLD:
      bond_broken = angle_energy(p1, p2, p3, iaparams, &ret);
      break; 
#endif
#ifdef BOND_ANGLE
    case BONDED_IA_ANGLE_HARMONIC:
      bond_broken = angle_harmonic_energy(p1, p2, p3, iaparams, &ret);
      break;
    case BONDED_IA_ANGLE_COSINE:
      bond_broken = angle_cosine_energy(p1, p2, p3, iaparams, &ret);
      break;
    case BONDED_IA_ANGLE_COSSQUARE:
      bond_broken = angle_cossquare_energy(p1, p2, p3, iaparams, &ret);
      break;
#endif
#ifdef BOND_ANGLEDIST
    case BONDED_IA_ANGLEDIST:
      bond_broken = angledist_energy(p1, p2, p3, iaparams, &ret);
      break;
#endif
#ifdef BOND_ENDANGLEDIST
    case BONDED_IA_ENDANGLEDIST:
      bond_broken = endangledist_pair_energy(p1, p2, iaparams, dx, &ret);
      break;
#endif
    case BONDED_IA_DIHEDRAL:
      bond_broken = dihedral_energy(p2, p1, p3, p4, iaparams, &ret);
      break;
#ifdef BOND_CONSTRAINT
    case BONDED_IA_RIGID_BOND:
      bond_broken = 0;
      ret = 0;
      break;
#endif
#ifdef TABULATED
    case BONDED_IA_TABULATED:
      switch(iaparams->p.tab.type) {
      case TAB_BOND_LENGTH:
	bond_broken = tab_bond_energy(p1, p2, iaparams, dx, &ret);
	break;
      case TAB_BOND_ANGLE:
	bond_broken = tab_angle_energy(p1, p2, p3, iaparams, &ret);
	break;
      case TAB_BOND_DIHEDRAL:
	bond_broken = tab_dihedral_energy(p2, p1, p3, p4, iaparams, &ret);
	break;
      default :
	errtxt = runtime_error(128 + ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{072 add_bonded_energy: tabulated bond type of atom %d unknown\n", p1->p.identity);
	return;
      }
      break;
#endif
#ifdef OVERLAPPED
    case BONDED_IA_OVERLAPPED:
      switch(iaparams->p.overlap.type) {
      case OVERLAP_BOND_LENGTH:
        bond_broken = overlap_bond_energy(p1, p2, iaparams, dx, &ret);
        break;
      case OVERLAP_BOND_ANGLE:
        bond_broken = overlap_angle_energy(p1, p2, p3, iaparams, &ret);
        break;
      case OVERLAP_BOND_DIHEDRAL:
        bond_broken = overlap_dihedral_energy(p2, p1, p3, p4, iaparams, &ret);
        break;
      default :
        errtxt = runtime_error(128 + ES_INTEGER_SPACE);
        ERROR_SPRINTF(errtxt,"{072 add_bonded_energy: overlapped bond type of atom %d unknown\n", p1->p.identity);
        return;
      }
      break;
#endif
#ifdef BOND_VIRTUAL
    case BONDED_IA_VIRTUAL_BOND:
      bond_broken = 0;
      ret = 0;
      break;
#endif
    default :
      errtxt = runtime_error(128 + ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{073 add_bonded_energy: bond type of atom %d unknown\n", p1->p.identity);
      return;
    }

    if (bond_broken) {
      switch (n_partners) {
      case 1: {
	char *errtext = runtime_error(128 + 2*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{083 bond broken between particles %d and %d} ",
		      p1->p.identity, p2->p.identity); 
	break;
      }
      case 2: {
	char *errtext = runtime_error(128 + 3*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{084 bond broken between particles %d, %d and %d} ",
		      p1->p.identity, p2->p.identity, p3->p.identity); 
	break;
      }
      case 3: {
	char *errtext = runtime_error(128 + 4*ES_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{085 bond broken between particles %d, %d, %d and %d} ",
		      p1->p.identity, p2->p.identity, p3->p.identity, p4->p.identity); 
	break;
      }
      }
      // bond broken, don't add whatever we find in the energy
      continue;
    }

    *obsstat_bonded(&energy, type_num) += ret;
  }
}

/** Calculate kinetic energies for one particle.
    @param p1 particle for which to calculate energies
*/
inline void add_kinetic_energy(Particle *p1)
{
#ifdef VIRTUAL_SITES
  if (ifParticleIsVirtual(p1)) return;
#endif

  /* kinetic energy */
  energy.data.e[0] += (SQR(p1->m.v[0]) + SQR(p1->m.v[1]) + SQR(p1->m.v[2]))*PMASS(*p1);

#ifdef ROTATION
#ifdef ROTATION_PER_PARTICLE
if (p1->p.rotation)
#endif
{
#ifdef ROTATIONAL_INERTIA
  /* the rotational part is added to the total kinetic energy;
     Here we use the rotational inertia  */

  energy.data.e[0] += (SQR(p1->m.omega[0])*p1->p.rinertia[0] +
		       SQR(p1->m.omega[1])*p1->p.rinertia[1] +
		       SQR(p1->m.omega[2])*p1->p.rinertia[2])*time_step*time_step;
#else
  /* the rotational part is added to the total kinetic energy;
     at the moment, we assume unit inertia tensor I=(1,1,1)  */
  energy.data.e[0] += (SQR(p1->m.omega[0]) + SQR(p1->m.omega[1]) + SQR(p1->m.omega[2]))*time_step*time_step;
#endif
}
#endif	
}

/*@}*/

#endif
