#ifndef STATISTICS_H
#define STATISTICS_H
/** \file statistics.h
    This file contains the code for simply statistics on the data.

    <b>Responsible:</b>
    <a href="mailto:mann@mpip-mainz.mpg.de">BAM</a>

*/

#include <tcl.h>
#include "particle_data.h"

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Particles' initial positions (needed for g1(t), g2(t), g3(t) in \ref analyze) */
extern float *partCoord_g, *partCM_g;

/** Particles' current configuration (updated if NULL, set to NULL by on_particle_change and on_integration_start) */
extern Particle *partCfg;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Implements the Tcl command mindist [<posx> <posy> <posz>] or mindist [<part_id> <r_catch>]. 
    Without any parameters it returns the minimal distance of two particles (in minimum image convention). 
    If the coordinates mindist <posx> <posy> <posz> are given, it returns the minimum distance of all particles to that position.
    If <part_id> <r_catch> is given, it returns the identities of all particles which are less than <r_catch> away from <part_id>.
*/
int mindist(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** Implements the Tcl command analyze <what> <N_P> <MPC> [<N_CI> [<N_pS> <N_nS>]] [-g] for basic analysis.
    <li> <what> == 0  returns the quadratic end-to-end-distance averaged over all polymers
    <li> <what> == 1  returns the radius of gyration averaged over all chains
    <li> <what> == 2  returns the hydrodynamic radius 
    <li> <what> == 3  returns the mean-square displacement g1(t) of a monomer,
                              the mean-square displacement g2(t) of in the center of gravity of the chain itself,
                              the motion of the center of mass g3(t)
                      as a tcl-list {g1(t) g2(t) g3(t)} on its first call or if invoked with option -g it writes
                      the current configuration to float partCoord_g[3*n_total_particles] and float partCM_g[3*N_P] */
int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv);


#endif
