/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef STATISTICS_MOLECULE_H
#define STATISTICS_MOLECULE_H
/** \file statistics_molecule.hpp

    This file contains the code for statistics on the data using the
    molecule information set with analyse set, as it is described in
    the file \ref topology.hpp.
*/

#include "utils.hpp"
#include "statistics.hpp"
#include "topology.hpp"

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

#endif
