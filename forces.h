/**************************************************/
/*******************  FORCES.H  *******************/
/**************************************************/
#ifndef FORCES_H
#define FORCES_H
/** \file forces.h
    Force calculation. Header file for \ref forces.c "forces.c".
 */
/************************************************
 * exported variables
 ************************************************/

/** \name Exported Variables */
/*@{*/
/** the Bjerrum length for electrostatics */
extern double Bjerrum;
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
/*@}*/

#endif
