// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file energy.h
    Implementation of the energy calculation.
*/
#ifndef ENERGY_H
#define ENERGY_H
#include "utils.h"
#include "integrate.h"
#include "statistics.h"
#include "thermostat.h"
#include "communication.h"
#ifdef NPT
#include "pressure.h"
#endif

/* include the energy files */
#include "p3m.h"
#include "lj.h"
#include "buckingham.h"
#include "soft_sphere.h"
#include "ljcos.h"
#include "tab.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "subt_lj.h"
#include "angle.h"
#include "dihedral.h"
#include "debye_hueckel.h"
#include "mmm1d.h"
#include "mmm2d.h"
#include "maggs.h"
#include "morse.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat energy, total_energy;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** parallel energy calculation.
    @param result non-zero only on master node; will contain the cumulative over all nodes. */
void energy_calc(double *result);

/** Calculate non bonded energies between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param d         vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
MDINLINE void add_non_bonded_pair_energy(Particle *p1, Particle *p2, double d[3],
					 double dist, double dist2)
{
  IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);
  double ret = 0;

#ifdef LENNARD_JONES
  /* lennard jones */
  ret += lj_pair_energy(p1,p2,ia_params,d,dist);
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

#ifdef TABULATED
  /* tabulated */
  ret += tabulated_pair_energy(p1,p2,ia_params,d,dist);
#endif

#ifdef LJCOS
  /* lennard jones cosine */
  ret += ljcos_pair_energy(p1,p2,ia_params,d,dist);
#endif
  
#ifdef ROTATION
  /* Gay-Berne */
  ret += gb_pair_energy(p1,p2,ia_params,d,dist);
#endif
  *obsstat_nonbonded(&energy, p1->p.type, p2->p.type) += ret;

#ifdef ELECTROSTATICS
  if (coulomb.method != COULOMB_NONE) {
    /* real space coulomb */
    switch (coulomb.method) {
    case COULOMB_P3M:
      ret = p3m_coulomb_pair_energy(p1,p2,d,dist2,dist);
#ifdef DIPOLES
      ret += p3m_dipol_pair_energy(p1,p2,d,dist2,dist); 
#endif
      break;
    case COULOMB_DH:
      ret = dh_coulomb_pair_energy(p1,p2,dist);
      break;
    case COULOMB_MMM1D:
      ret = mmm1d_coulomb_pair_energy(p1,p2,d,dist2,dist);
      break;
    case COULOMB_MMM2D:
      ret = mmm2d_coulomb_pair_energy(p1,p2,d,dist2,dist);
      break;
    default :
      ret = 0.;
    }
    energy.coulomb[0] += ret;
  }
#endif
}

/** Calculate bonded energies for one particle.
    @param p1 particle for which to calculate energies
*/
MDINLINE void add_bonded_energy(Particle *p1)
{
  char *errtxt;
  Particle *p2, *p3 = NULL, *p4 = NULL;
  Bonded_ia_parameters *iaparams;
  int i, type_num, type, n_partners, bond_broken;
  double ret, dx[3] = {0, 0, 0};

  i = 0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i++];
    iaparams = &bonded_ia_params[type_num];
    type = iaparams->type;
    n_partners = iaparams->num;
    
    /* fetch particle 2, which is always needed */
    p2 = local_particles[p1->bl.e[i++]];
    if (!p2) {
      errtxt = runtime_error(128 + 2*TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{069 bond broken between particles %d and %d (particles not stored on the same node)} ",
	      p1->p.identity, p1->bl.e[i-1]);
      return;
    }

    /* fetch particle 3 eventually */
    if (n_partners >= 2) {
      p3 = local_particles[p1->bl.e[i++]];
      if (!p3) {
	errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{070 bond broken between particles %d, %d and %d (particles not stored on the same node)} ",
		p1->p.identity, p1->bl.e[i-2], p1->bl.e[i-1]);
	return;
      }
    }

    /* fetch particle 4 eventually */
    if (n_partners >= 3) {
      p4 = local_particles[p1->bl.e[i++]];
      if (!p4) {
	errtxt = runtime_error(128 + 4*TCL_INTEGER_SPACE);
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
    case BONDED_IA_ANGLE:
      bond_broken = angle_energy(p1, p2, p3, iaparams, &ret);
      break;
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
	errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtxt,"{072 add_bonded_energy: tabulated bond type of atom %d unknown\n", p1->p.identity);
	return;
      }
      break;
#endif
    default :
      errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"{073 add_bonded_energy: bond type of atom %d unknown\n", p1->p.identity);
      return;
    }

    switch (n_partners) {
    case 1:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 2*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{074 bond broken between particles %d and %d} ",
		p1->p.identity, p2->p.identity);
	continue;
      }
      break;
    case 2:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 3*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{075 bond broken between particles %d, %d and %d} ",
		p1->p.identity, p2->p.identity, p3->p.identity); 
	continue;
      }
      break;
    case 3:
      if (bond_broken) {
	char *errtext = runtime_error(128 + 4*TCL_INTEGER_SPACE);
	ERROR_SPRINTF(errtext,"{076 bond broken between particles %d, %d, %d and %d} ",
		p1->p.identity, p2->p.identity, p3->p.identity, p4->p.identity); 
	continue;
      }
    }

    *obsstat_bonded(&energy, type_num) += ret;
  }
}

/** Calculate kinetic energies for one particle.
    @param p1 particle for which to calculate energies
*/
MDINLINE void add_kinetic_energy(Particle *p1)
{
  /* kinetic energy */
  energy.data.e[0] += (SQR(p1->m.v[0]) + SQR(p1->m.v[1]) + SQR(p1->m.v[2]))*PMASS(*p1);

#ifdef ROTATION
  /* the rotational part is added to the total kinetic energy;
     at the moment, we assume unit inertia tensor I=(1,1,1)  */
  energy.data.e[0] += (SQR(p1->m.omega[0]) + SQR(p1->m.omega[1]) + SQR(p1->m.omega[2]))*time_step*time_step;
#endif	
}

/** implementation of analyze energy */
int parse_and_print_energy(Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
