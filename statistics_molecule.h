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

/** Using the topology information stored in \ref topology::topology
    this routine foldes all particles belonging to a molecule such
    that the center of mass of the molecule is inside the simulation
    box.  

@param coord is an array specifying the full coordinates of all
particles.  Particles must be stored in particle ID order in this
array. @param shift is a vector specifying a shift for the entire
system coordinates that is applied prior to folding: This is used for
visualization purposes.

*/
int analyze_fold_molecules(float *coord, double shift[3]);

/* calculate the center of mass of a molecule */
void calc_mol_center_of_mass(Molecule mol, double com[3]);

/* calculate the center of mass of a molecule as above but uses mass of the particles*/
void mol_center_of_mass_(Molecule mol, double com[3]);


