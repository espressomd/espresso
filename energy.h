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
