// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef FORCES_H
#define FORCES_H
/** \file forces.h Force calculation. 
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  \todo Preprocessor switches for all forces (Default: everything is turned on).
 *  \todo Implement more flexible thermostat, e.g. which thermostat to use.
 *
 *  For more information see \ref forces.c "forces.c".
 */
#include <tcl.h>
#include "config.h"
#include "thermostat.h"

/* include the force files */
#include "p3m.h"
#include "lj.h"
#include "ljcos.h"
#include "gb.h"
#include "fene.h"
#include "harmonic.h"
#include "angle.h"
#include "debye_hueckel.h"
#include "mmm1d.h"
#include "constraint.h"

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** initialize real particle forces with thermostat forces and
    ghost particle forces with zero. */
void init_forces();

/** Calculate forces.
 *
 *  A short list, what the function is doing:
 *  <ol>
 *  <li> Initialize forces with: \ref friction_thermo (ghost forces with zero).
 *  <li> Calculate \ref tcl_bonded "bonded interaction" forces:<br>
 *       Loop all local particles (not the ghosts). 
 *       <ul>
 *       <li> FENE
 *       <li> ANGLE (cos bend potential)
 *       </ul>
 *  <li> Calculate \ref tcl_non_bonded "non-bonded short range interaction" forces:<br>
 *       Loop all \ref IA_Neighbor::vList "verlet lists" of all \ref #cells.
 *       <ul>
 *       <li> Lennard-Jones.
 *       <li> Real space part: Coulomb.
 *       <li> Ramp.
 *       </ul>
 *  <li> Calculate long range interaction forces:<br>
         Uses \ref P3M_calc_kspace_forces.
 *  </ol>
 */
void force_calc();

/** Calculate long range forces (P3M, MMM1D, MMM2d...). */
void calc_long_range_forces();

/** Calculate non bonded forces between a pair of particles.
    @param p1        pointer to particle 1.
    @param p2        pointer to particle 2.
    @param ia_params interaction parameters between p1 and p2. 
    @param d[3]      vector between p1 and p2. 
    @param dist      distance between p1 and p2.
    @param dist2     distance squared between p1 and p2. */
MDINLINE void add_non_bonded_pair_force(Particle *p1, Particle *p2, 
					double d[3], double dist, double dist2)
{
  IA_parameters *ia_params = get_ia_param(p1->p.type,p2->p.type);

  /* lennard jones */
  add_lj_pair_force(p1,p2,ia_params,d,dist);
  /* lennard jones cosine */
  add_ljcos_pair_force(p1,p2,ia_params,d,dist);
  
#ifdef ROTATION
  /* Gay-Berne */
  add_gb_pair_force(p1,p2,ia_params,d,dist);
#endif

#ifdef ELECTROSTATICS
  /* real space coulomb */
  if(coulomb.method == COULOMB_P3M) 
    add_p3m_coulomb_pair_force(p1,p2,d,dist2,dist);
  else if(coulomb.method == COULOMB_DH)
    add_dh_coulomb_pair_force(p1,p2,d,dist);
#endif

}

/** Calculate bonded forces between a pair of particles.
    Also contains the constraint forces, which are treated
    as bonded interactions in Espresso.
    @param p1 particle for which to calculate forces
*/
MDINLINE void add_bonded_pair_force(Particle *p1)
{
  int i, type_num;

  i=0;
  while(i<p1->bl.n) {
    type_num = p1->bl.e[i];
    switch(bonded_ia_params[type_num].type) {
    case BONDED_IA_FENE:
      add_fene_pair_force(p1,
			  checked_particle_ptr(p1->bl.e[i+1]), type_num);
      i+=2; break;
    case BONDED_IA_HARMONIC:
      add_harmonic_pair_force(p1,
			      checked_particle_ptr(p1->bl.e[i+1]), type_num);
      i+=2; break;
    case BONDED_IA_ANGLE:
      add_angle_force(p1,
		      checked_particle_ptr(p1->bl.e[i+1]),
		      checked_particle_ptr(p1->bl.e[i+2]), type_num);
      i+=3; break;
    default :
      fprintf(stderr,"WARNING: Bonds of atom %d unknown\n",p1->p.identity);
      i = p1->bl.n; 
      break;
    }
  }
  
#ifdef CONSTRAINTS
  add_constraints_forces(p1);
#endif
}

/** add force to another. This is used when collecting ghost forces. */
MDINLINE void add_force(ParticleForce *F_to, ParticleForce *F_add)
{
  int i;
  for (i = 0; i < 3; i++)
    F_to->f[i] += F_add->f[i];
#ifdef ROTATION
  for (i = 0; i < 3; i++)
    F_to->torque[i] += F_add->torque[i];
#endif
}

/*@}*/

#endif
