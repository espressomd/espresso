/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you
 * accepted upon receiving the distribution and by which you are
 * legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * You should have received a copy of that license along with this
 * program; if not, refer to http://www.espresso.mpg.de/license.html
 * where its current version can be found, or write to
 * Max-Planck-Institute for Polymer Research, Theory Group, 
 * PO Box 3148, 55021 Mainz, Germany. 
 * Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
 */

/** \file statistics_fluid.h
 *
 * Fluid related analysis functions.
 * Header file for \ref statistics_fluid.c.
 *
 */

#ifndef STATISTICS_FLUID_H
#define STATISTICS_FLUID_H

#include "utils.h"

#ifdef LB

/** Caclulate mass of the LB fluid.
 * \param result Fluid mass
 */
void lb_calc_fluid_mass(double *result);

/** Calculate momentum of the LB fluid.
 * \param result Fluid momentum
 */
void lb_calc_fluid_momentum(double *result);

/** Calculate temperature of the LB fluid.
 * \param result Fluid temperature
 */
void lb_calc_fluid_temp(double *result);

/** Parser for fluid related analysis functions. */
int parse_analyze_fluid(Tcl_Interp *interp, int argc, char **argv);

#endif /* LB */

#endif /* STATISTICS_FLUID_H */
