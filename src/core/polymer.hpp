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
    configuration of (partially charged) polymer chains with
    counterions and salt molecules, assigning velocities to the
    particles and cross-linking the polymers if necessary.

    For more information on polymer, see polymer.cpp.
*/

#include "PartCfg.hpp"
#include "particle_data.hpp"
#include "utils/Vector.hpp"

Utils::Vector3d random_position(const std::function<double()> &generate_rn);
Utils::Vector3d random_unit_vector(const std::function<double()> &generate_rn);

/** Returns the miminum distance between position pos and all existing
 * particles.
 */
double mindist(PartCfg &partCfg, const Utils::Vector3d pos);

/** Determines whether a given position pos is valid, i.e., it doesn't collide
 * with existing or buffered particles, nor with existing constraints (if
 * respect_constraints).
 */
bool is_valid_position(
    const Utils::Vector3d *pos,
    const std::vector<std::vector<Utils::Vector3d>> *positions,
    const PartCfg &partCfg, const double min_distance,
    const int respect_constraints);

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
draw_polymer_positions(PartCfg &partCfg, const int n_polymers,
                       const int beads_per_chain, double const bond_length,
                       const std::vector<Utils::Vector3d> &start_positions,
                       const double min_distance, const int max_tries,
                       const int use_bond_angle, const double bond_angle,
                       const int respect_constraints, const int seed);

#endif
