#ifndef STATISTICS_H
#define STATISTICS_H
/** \file statistics.h
    This file contains the code for simply statistics on the data.

    <b>Responsible:</b>
    <a href="mailto:svanebor@mpip-mainz.mpg.de">Carsten</a>

*/

#include <tcl.h>

/** Implements the Tcl command mindist. It returns the minimal distance of two particles
    (in minimum image convention). */
int mindist(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implements the Tcl command analyze <what> <N_P> <MPC> [<N_CI> [<N_pS> <N_nS>]] for basic analysis.
    <li> <what> == 0  returns the quadratic end-to-end-distance averaged over all polymers
    <li> <what> == 1  returns the radius of gyration averaged over all chains
    <li> <what> == 2  returns the hydrodynamic radius */
int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv);


#endif
