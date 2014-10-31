/*
  Copyright (C) 2014 The ESPResSo project
  
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
/* vim: set ts=8 sts=2 sw=2 et: */

#ifndef __READ_PDB_HPP
#define __READ_PDB_HPP

#include "PdbParser.hpp"
#include "particle_data.hpp"

/** Call only on the master node: Parse pdb file and add contained particles.
    @param pdb_file Filename of the pdb file.
    @first_id Id of the first particle to add.
    @type Type for the particles.
    @return Number of particles that were added.
 */

int pdb_add_particles_from_file(char *pdb_file, int first_id, int type);
#endif
