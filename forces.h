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

/* include the force files */
#include "p3m.h"
#include "lj.h"
#include "ljcos.h"
#include "fene.h"
#include "harmonic.h"
#include "angle.h"
#include "debye_hueckel.h"
#include "mmm1d.h"
#include "constraint.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** The minimum particle distance seen by the ramp potential. */
extern double minimum_part_dist;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Calculates lj cap radii */
void force_init();

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

/** Calculate bonded interactions (forces) for Particle List p (length np).
    This includes also the constraint forces.
    @param p Particle List
    @param np length of that list
*/
void calc_bonded_forces(Particle *p, int np);

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
					IA_parameters *ia_params, 
					double d[3], double dist, double dist2)
{
  /* lennard jones */
  add_lj_pair_force(p1,p2,ia_params,d,dist);
  /* lennard jones cosine */
  add_ljcos_pair_force(p1,p2,ia_params,d,dist);

#ifdef ELECTROSTATICS
  /* real space coulomb */
  if(coulomb.method == COULOMB_P3M) 
    add_p3m_coulomb_pair_force(p1,p2,d,dist2,dist);
  else if(coulomb.method == COULOMB_DH)
    add_dh_coulomb_pair_force(p1,p2,d,dist);
#endif

}

/*@}*/

#endif
