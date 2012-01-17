/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef EWALD_H 
#define EWALD_H
/** \file ewald.h   Ewald algorithm for long range coulomb interaction.
 *
 *  Implementation of the standard Ewald-Summation
 *
 *  Further reading: 
 *  <ul>
 *  <li> P.P. Ewald,
 *       <i>Die Berechnung optischer und elektrostatischer Gitterpotentiale</i>,
 *       Ann. Phys. (64) 253-287, 1921
 *  <li> M. Deserno and C. Holm,
 *       <i>How to mesh up {E}wald sums. I. + II.</i>,
 *       J. Chem. Phys. (109) 7678, 1998; (109) 7694, 1998
 *  <li> M. Deserno, C. Holm and H. J. Limbach,
 *       <i>How to mesh up {E}wald sums. </i>,
 *       in Molecular Dynamics on Parallel Computers,
 *       Ed. R. Esser et al., World Scientific, Singapore, 2000
 *  <li> 
 *  </ul>
 *
 *  For more information about the ewald algorithm,
 *  see \ref ewald.c "ewald.c"
 */

#include "config.h"
#include "debug.h"

#include "interaction_data.h"

#ifdef ELECTROSTATICS

/** This value for ewald.epsilon indicates metallic boundary conditions. */
#define EWALD_EPSILON_METALLIC 0.0

/************************************************
 * data types
 ************************************************/

/** Structure to hold Ewald parameters and some dependend variables. */
typedef struct {
  /** Ewald splitting parameter (0<alpha<1), rescaled to alpha_L = alpha * box_l. */
  double alpha_L;
  /** Cutoff radius for real space electrostatics (>0), rescaled to r_cut_iL = r_cut * box_l_i. */
  double r_cut_iL;
  /** epsilon of the "surrounding dielectric". */
  double epsilon;
  /** unscaled \ref alpha_L for use with fast inline functions only */
  double alpha;
  /** unscaled \ref r_cut_iL for use with fast inline functions only */
  double r_cut;
  /** Maximum KVec, e.g. maximum frequency in k-space*/
  int kmax;
  /** squared \ref kmax */
  int kmaxsq;
} ewald_struct;

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Ewald parameters. */
extern ewald_struct ewald;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/// print the ewald parameters to the interpreters result
int tclprint_to_result_EWALD(Tcl_Interp *interp);

/// parse the ewald parameters
int tclcommand_inter_coulomb_parse_ewald(Tcl_Interp * interp, int argc, char ** argv);
  
/// sanity checks
int EWALD_sanity_checks();

/** Initialize all structures, parameters and arrays needed for the 
 *  Ewald algorithm.
 */
void EWALD_init();

/** Calculate number of charged particles, the sum of the squared
    charges and the squared sum of the charges. */
void EWALD_count_charged_particles();

/** Reallocate memory for k-space caches */
void EWALD_on_resort_particles();

/** Updates \ref ewald_struct::alpha and \ref ewald_struct::r_cut if \ref box_l changed. */
void EWALD_scaleby_box_l();

/** Calculate the k-space contribution to the coulomb interaction forces. */ 
double EWALD_calc_kspace_forces(int force_flag, int energy_flag);

/** Calculate real space contribution of coulomb pair forces.
    Returns the energy, which is needed for NPT. */
MDINLINE double add_ewald_coulomb_pair_force(Particle *p1, Particle *p2,
					   double *d,double dist2,double dist,double force[3])
{
  int j;
  double fac1,fac2, adist, erfc_part_ri;
#ifdef ONEPART_DEBUG
  double fac=0.; /*TODO: (BUG) fac was not declared, now code compiles but its value (0) is arbitrary */
#endif 
  if(dist < ewald.r_cut) {
    adist = ewald.alpha * dist;
#if USE_ERFC_APPROXIMATION
    erfc_part_ri = AS_erfc_part(adist) / dist;
    fac1 = coulomb.prefactor * p1->p.q * p2->p.q  * exp(-adist*adist);
    fac2 = fac1 * (erfc_part_ri + 2.0*ewald.alpha*wupii) / dist2;
#else
    erfc_part_ri = erfc(adist) / dist;
    fac1 = coulomb.prefactor * p1->p.q * p2->p.q;
    fac2 = fac1 * (erfc_part_ri + 2.0*ewald.alpha*wupii*exp(-adist*adist)) / dist2;
#endif
    for(j=0;j<3;j++)
      force[j] += fac2 * d[j];
    ESR_TRACE(fprintf(stderr,"%d: RSE: Pair (%d-%d) dist=%.3f: force (%.3e,%.3e,%.3e)\n",this_node,
 		      p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: ESR  f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac2));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: ESR  f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac2));

    return fac1 * erfc_part_ri;
  }
  return 0.0;
}

/** Calculate real space contribution of coulomb pair energy. */
MDINLINE double ewald_coulomb_pair_energy(Particle *p1, Particle *p2,
				     double *d,double dist2,double dist)
{
  double adist, erfc_part_ri;

  if(dist < ewald.r_cut) {
    adist = ewald.alpha * dist;
#if USE_ERFC_APPROXIMATION
    erfc_part_ri = AS_erfc_part(adist) / dist;
    return coulomb.prefactor*p1->p.q*p2->p.q *erfc_part_ri*exp(-adist*adist);
#else
    erfc_part_ri = erfc(adist) / dist;
    return coulomb.prefactor*p1->p.q*p2->p.q *erfc_part_ri;
#endif
  }
  return 0.0;
}

/** Clean up Ewald memory allocations. */
void   Ewald_exit();

/*@}*/
#endif

#endif
