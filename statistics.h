#ifndef STATISTICS_H
#define STATISTICS_H
/** \file statistics.h
    This file contains the code for simply statistics on the data.

    <b>Responsible:</b>
    <a href="mailto:mann@mpip-mainz.mpg.de">BAM</a>

*/

#include <tcl.h>
#include "particle_data.h"

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Implements the Tcl command analyze <what> <structure info> ... for basic analysis.
    Possible structures are
    <ul>
    <li> chain a set of equal-length chains. Has structure info <chain_start> <n_chains> <chain_length>
    </ul>
    <li> <what> == set   define the structure. The second argument defines the topology
    to set, i. e. chain at the moment.
    <li> <what> == mindist returns the minimal distance of two particles (needs no structure info)
    <li> <what> == nbhood returns all particles within a given radius around a fixed particle.
    (needs the fixed particles id
    and the radius of the neighborhood)
    <li> <what> == distto returns the minimal distance of a particle to a fixed point
    (needs that points coordinates)
    <li> <what> == re    returns the quadratic end-to-end-distance averaged over all polymers
    (chain structure)
    <li> <what> == rg    returns the radius of gyration averaged over all chains (chain structure)
    <li> <what> == rh    returns the hydrodynamic radius  (chain structure)
    <li> <what> == g123  returns the mean-square displacement g1(t) of a monomer,
                              the mean-square displacement g2(t) of in the center
			      of gravity of the chain itself,
                              the motion of the center of mass g3(t)
			 as a tcl-list {g1(t) g2(t) g3(t)} (needs chain structure).
			 If before the structure info you give -init,
			 the current configuration is stored as reference config
    All tasks below the end-to-end distance need the particles to be stored
    consecutively starting with identity 0.
*/
int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/*@}*/

#endif
