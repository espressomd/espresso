// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef STATISTICS_MOLECULE_H
#define STATISTICS_MOLECULE_H
/** \file statistics_molecule.h 

    This file contains the code for statistics on the data using the
    molecule information set with analyse set, as it is described in 
    the file \ref topology.h.

    <b>Responsible:</b>
    <a href="mailto:hanjo@mpip-mainz.mpg.de">Hanjo</a>
*/

#endif
#include "statistics.h"
#include "parser.h"
#include "debug.h"
#include "topology.h"

/** Using the topology information stored in \ref molecules this routine foldes all particles belonging to a molecule 
 Using the values of \ref chain_start , \ref chain_length and \ref
chain_n_chains this routine folds a set of chains of uniform length
without losing bonding connectivity.  For this to work the chains must
be of uniform length and all the bonds between two particles must be
specified as belonging to the particle with the lower particle ID. \n
For example a bond should be assigned in this way: \n
<tt> part \<partnum\> bond \<bond_type\> \<partnum + 1\> </tt> \n

All bonded particles must also have sequential identities corresponding to sequence along the chain.

@param coord is an array specifying the full coordinates of all particles.

*/
int analyze_fold_chains(float *coord);

/* calculate the center of mass of a molecule */
void calc_center_of_mass(Molecule mol, double com[3]);
