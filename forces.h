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

/*@}*/

#endif
