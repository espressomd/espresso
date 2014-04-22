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

#include <mpi.h>
#include "utils.hpp"
#include "grid.hpp"

/** \file molforces.hpp
 *  Routines for calculating and applying trap forces upon molecules.
 *  This trap force can be set to
 *  - a harmonic potential with a restlength of zero on the molecular centre of mass
 *  - a drag on the molecular velocity
 *  - a cancelation of the total force on the molecule (including thermostat forces)
 *  The centre of mass can be fixed to an absolute position or to a relative position in the
 *  simulation box.
 *  The molecular trap forces is distributed evenly upon all particles in a molecule.
 *  (see file \ref topology.cpp and file \ref molforces.cpp)  
 */

#ifdef MOLFORCES

extern int     IsTrapped;

/**
   Checks if there are any molecules trapped (IsTrapped=1) and if so calls calc_mol_info
   and apply_mol_constaints */
void calc_and_apply_mol_constraints();

#endif
