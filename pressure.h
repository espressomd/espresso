#ifndef PRESSURE_H
#define PRESSURE_H
#include "statistics.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/
///
extern Observable_stat virials;
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initializes extern Energy_stat \ref #virials to be used by \ref calc_virials. */
void init_virials();

/** Calculates the virials of the system in parallel (hence it should be called by \ref mpi_gather_stats with job=2).<BR>
    Due to the nature of a virial being <tt>Sum(i=0..n_total_particles)(Sum(j=i+1..n_total_particles)(r_ij*F_ij))</tt>
    this function is based on a merge of \ref force_calc into \ref calc_energy. */
void calc_virials();

/** Calculates the pressure in the system from a virial expansion using the terms from \ref calc_virials.<BR>
    Output is stored in the \ref #virials array, in which (on the first node) each component carries the corresponding pressure,
    while <tt>virials.sum.e[0]</tt> contains the total pressure, <tt>virials.node.e[0]</tt> the sum of all squared pressure components,
    <tt>virials.sum.e[1]</tt> the pressure of the ideal gas, and <tt>virials.node.e[1]</tt> the kinetic energy.
*/
void calc_pressure(void);

/** implementation of analyze pressure */
int parse_and_print_pressure(Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
