/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef CORE_TABULATED_HPP
#define CORE_TABULATED_HPP

/** \file nonbonded_tab.hpp
 *  Routines to calculate the  energy and/or  force
 *  for a particle pair or bonds via interpolating from lookup tables.
 *  \ref forces.cpp
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 */

#include "config.hpp"

#ifdef TABULATED

#include "debug.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

/** Non-Bonded tabulated potentials:
    Reads tabulated parameters and force and energy tables from a file.
    ia_params and force/energy tables are then communicated to each node

    @param part_type_a particle type for which the interaction is defined
    @param part_type_b particle type for which the interaction is defined
    @param filename from which file to fetch the data

    @return <ul>
    <li> 0 on success
    <li> 1 on particle type mismatches
    <li> 2 file name too long
    <li> 3 cannot open the file
    <li> 4 file too short
    <li> 5 file broken, cannot parse numbers
    <li> 6 number of points of existing potential changed
    </ul>
*/
int tabulated_set_params(int part_type_a, int part_type_b, double min,
                         double max, std::vector<double> const &energy,
                         std::vector<double> const &force);

/** Add a non-bonded pair force by linear interpolation from a table.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline void add_tabulated_pair_force(const Particle *const p1,
                                     const Particle *const p2,
                                     IA_parameters const *ia_params,
                                     double const d[3], double dist,
                                     double force[3]) {
  if (dist < ia_params->TAB.cutoff()) {
    auto const fac = ia_params->TAB.force(dist) / dist;

    for (int j = 0; j < 3; j++)
      force[j] -= fac * d[j];
  }
}

/** Add a non-bonded pair energy by linear interpolation from a table.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
inline double tabulated_pair_energy(Particle const *, Particle const *,
                                    IA_parameters const *ia_params,
                                    const double d[3], double dist) {
  if (dist < ia_params->TAB.cutoff()) {
    return ia_params->TAB.energy(dist);
  } else {
    return 0.0;
  }
}

#endif

#endif
