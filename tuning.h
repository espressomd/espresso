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
