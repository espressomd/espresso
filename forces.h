/**************************************************/
/*******************  FORCES.H  *******************/
/**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "communication.h"
#include "global.h"


/*******************  Functions  *******************/

/** initialization of force calculation. 
 *  init interaction matrices for short range forces. 
 *  init long range forces (P3M, MMM, etc).
*/
void force_init();

/** Calculate forces (and energies). */
void force_calc();

/** clean up the force part.*/
void force_exit(); 

