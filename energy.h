// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef ENERGY_H
#define ENERGY_H
#include "integrate.h"
#include "config.h"
#include "statistics.h"
#include "thermostat.h"
#include "communication.h"
#ifdef NPT
#include "pressure.h"
#endif

/* include the energy files */
#include "p3m.h"
#include "lj.h"
#include "ljcos.h"
#include "tab.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "subt_lj_harm.h"
#include "subt_lj_fene.h"
#include "angle.h"
#include "debye_hueckel.h"
#include "mmm1d.h"
#include "mmm2d.h"


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
  ret  += lj_pair_energy(p1,p2,ia_params,d,dist);
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
  if (coulomb.bjerrum != 0.0) {
    /* real space coulomb */
    switch (coulomb.method) {
    case COULOMB_P3M:
      ret = p3m_coulomb_pair_energy(p1,p2,d,dist2,dist);
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
  int i, type_num;
  double ret;
  i=0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i];
    switch(bonded_ia_params[type_num].type) {
    case BONDED_IA_FENE:
      ret = fene_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
      i+=2; break;
    case BONDED_IA_HARMONIC:
      ret = harmonic_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
      i+=2; break;
    case BONDED_IA_SUBT_LJ_HARM:
      ret = subt_lj_harm_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
      i+=2; break; 
    case BONDED_IA_SUBT_LJ_FENE:
      ret = subt_lj_fene_pair_energy(p1, checked_particle_ptr(p1->bl.e[i+1]), type_num);
      i+=2; break; 
    case BONDED_IA_ANGLE:
      ret = angle_energy(p1,
			 checked_particle_ptr(p1->bl.e[i+1]),
			 checked_particle_ptr(p1->bl.e[i+2]), type_num);
      i+=3; break;
    default :
      fprintf(stderr,"add_bonded_energy: WARNING: Bonds of atom %d unknown\n",p1->p.identity);
      ret = 0;
      i = p1->bl.n;
      break;
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
  energy.data.e[0] += SQR(p1->m.v[0]) + SQR(p1->m.v[1]) + SQR(p1->m.v[2]);

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
