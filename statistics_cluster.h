// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef STATISTICS_CLUSTER_H
#define STATISTICS_CLUSTER_H
/** \file statistics_cluster.h
 *
 *  This file contains the necklace cluster algorithm. It can be used
 *  to identify the substructures 'pearls' and 'strings' on a linear
 *  chain.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 */

#include <tcl.h>
#include "particle_data.h"
#include "parser.h"

/** Parser for the necklace cluster algorithm

    \verbatim analyze necklace <pearl_treshold> <back_dist> <space_dist> <first> <length> \endverbatim 

    For more information see: Limbach H.J. and Holm C.  Single-Chain
    Properties of Polyelectrolytes in Poor Solvent J. Phys. Chem. B,
    107 (32), 8041 -8055, AUG 14 2003.  The first three parameters are
    tune parameters for the algorithm: pearl_treshold is the minimal
    number of monomers in a pearl. back_dist is the number of monomers
    along the chain backbone which are excluded from the space
    distance criterion to form clusters. space_dist is the distance
    between two monomers up to which they are considered to belong to
    the same clusters. The three parameters may be connected by
    scaling arguments. Make sure that your results are only weakly
    dependent on the exact choice of your parameters. For the
    algorithm the coordinates stored in \ref partCfg are used. The
    chain itself is defined by the identity first of its first monomer
    and the chain length length.
*/
int parse_necklace_analyzation(Tcl_Interp *interp, int argc, char **argv);
#endif
