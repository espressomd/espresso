#ifndef ENERGY_H
#define ENERGY_H
#include "statistics.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat energy;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

void init_energies();

void calc_energy();

/** implementation of analyze energy */
int parse_and_print_energy(Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
