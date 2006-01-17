// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file tuning.h
    This contains a timing loop for the force calculation. Via the global variable timings you can specify how many
    force evaluations are sampled. Via \ref markTime and \ref diffTime you can also easily time anything other than
    the force evaluation.
*/

#ifndef TUNING_H
#define TUNING_H

/** if positive, the number of samples for timing */
extern int timing_samples;

/** returns the time for some force calculations.
    Actually performs \ref mpi_integrate (0)
    @param default_samples the number of samples to take if
    \ref timing_samples is not set. */
double time_force_calc(int default_samples);

/** callback for \ref timing_samples */
int timings_callback(Tcl_Interp *interp, void *data);

/** set a time marker. \ref diffTime always gives the time in ms between
    the last two calls to markTime. */
void markTime();

/** calculate milliseconds between last two calls to \ref markTime. */
double diffTime();

#endif
