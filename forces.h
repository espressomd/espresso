/**************************************************/
/*******************  FORCES.H  *******************/
/**************************************************/
#ifndef FORCES_H
#define FORCES_H
/** \file forces.h
    Force calculation. Header file for \ref forces.c "forces.c".
 */
#include <tcl.h>

/************************************************
 * exported variables
 ************************************************/

/** \name Exported Variables */
/*@{*/
/** the Bjerrum length for electrostatics */
extern double Bjerrum;
/** the minimum particle distance seen by the ramp potential. */
extern double minimum_part_dist;
/*@}*/

/*******************  Functions  *******************/

/** \name Exported Functions */
/*@{*/

/** Initialization of force calculation. 
    init interaction matrices for short range forces. 
    init long range forces (P3M, MMM, etc).
*/
void force_init();
 

/** Calculate forces (and energies). */
void force_calc();

/** Clean up the force part. */
void force_exit(); 

/** Callback for setmd bjerrum. If the Bjerrum length is 0, the cutoff is also set to 0
    to avoid unnecessary computation of the electrostatics. */
int bjerrum_callback(Tcl_Interp *interp, void *_data);

/*@}*/

#endif
