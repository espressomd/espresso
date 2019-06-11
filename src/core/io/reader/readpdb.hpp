/*
  Copyright (C) 2014-2018 The ESPResSo project

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

#ifndef __READ_PDB_HPP
#define __READ_PDB_HPP

#include "PdbParser.hpp"
#include "particle_data.hpp"

namespace Reader {
namespace PDB {

struct PdbLJInteraction {
  int other_type;
  double epsilon, sigma;
};

/** Call only on the master node: Parse pdb file and add contained particles.
    @param pdb_file Filename of the pdb file.
    @param first_id Id of the first particle to add.
    @param type Type for the particles.
    @param ljInteractions LJ interactions vector.
    @param lj_rel_cutoff LJ cutoff, as a multiple of the read LJ sigma value.
    @param itp_file Optional ITP file.
    @param first_type Number of atom types to skip in the ITP file.
    @param fit Should the box be rescaled to hold the particles.
    @param lj_internal Should LJ interactions within the molecule be added.
    @param lj_diagonal Just the diagonal interaction terms of @p lj_internal.
    @return Number of particles that were added.
 */
int pdb_add_particles_from_file(char *pdb_file, int first_id, int type,
                                std::vector<PdbLJInteraction> &ljInteractions,
                                double lj_rel_cutoff = 2.5,
                                char *itp_file = nullptr, int first_type = 0,
                                bool fit = false, bool lj_internal = false,
                                bool lj_diagonal = false);
} // namespace PDB
} // namespace Reader
#endif
