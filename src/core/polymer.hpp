/*
  Copyright (C) 2010-2018 The ESPResSo project
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
#ifndef POLYMER_H
#define POLYMER_H
/** \file

    This file contains everything needed to create a start-up
    configuration of polymer chains which may respect already
    existing particles and/or constraints.

*/

#include "PartCfg.hpp"
#include "particle_data.hpp"
#include "utils/Vector.hpp"

Utils::Vector3d random_position(std::function<double()> const &generate_rn);
Utils::Vector3d random_unit_vector(std::function<double()> const &generate_rn);

/** Returns the miminum distance between position pos and all existing
 *  particles.
 */
double mindist(PartCfg &partCfg, Utils::Vector3d const &pos);

/** Determines whether a given position pos is valid, i.e., it doesn't collide
 *  with existing or buffered particles, nor with existing constraints (if
 *  respect_constraints). Returns true if valid, false if not.
 *  @param pos                 = the trial position in question <br>
 *  @param positions           = buffered positions to respect <br>
 *  @param partCfg             = existing particles to respect <br>
 *  @param min_distance        = threshold for the minimum distance between
 *  trial position and buffered/existing particles <br>
 *  @param respect_constraints = whether to respect constraints <br>
 */
bool is_valid_position(
    Utils::Vector3d const *pos,
    std::vector<std::vector<Utils::Vector3d>> const *positions,
    PartCfg const &partCfg, double min_distance, int respect_constraints);

/** Determines valid polymer positions and returns them.
 *  @param  n_polymers         = how many polymers to create <br>
 *  @param  beads_per_chain    = monomers per chain <br>
 *  @param  bond_length        = length of the bonds between two monomers <br>
 *  @param  seed               = seed for RNG <br>
 *  <br>
 *  @param  min_distance       = minimum distance between all particles <br>
 *  @param  max_try            = how often a monomer/polymer should be reset if
 *  current position collides with a previous particle <br>
 *  @param  angle              = desired bond-angle to be fixed <br>
 *  @param  respect_constrains = shall constraints be respected when setting up
 *  polymer?  (0=no, 1=yes, default: 0)
 */
std::vector<std::vector<Utils::Vector3d>>
draw_polymer_positions(PartCfg &partCfg, int n_polymers, int beads_per_chain,
                       double bond_length,
                       std::vector<Utils::Vector3d> const &start_positions,
                       double min_distance, int max_tries, int use_bond_angle,
                       double bond_angle, int respect_constraints, int seed);

#endif
