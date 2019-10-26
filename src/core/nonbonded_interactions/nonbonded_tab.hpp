/*
 * Copyright (C) 2010-2016,2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef CORE_TABULATED_HPP
#define CORE_TABULATED_HPP

/** \file
 *  Routines to calculate the energy and/or force for particle pairs via
 *  interpolation of lookup tables.
 *
 *  Implementation in \ref nonbonded_tab.cpp.
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 */

#include "config.hpp"

#ifdef TABULATED

#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

/** Set the parameters of a non-bonded tabulated potential.
 *  ia_params and force/energy tables are communicated to each node
 *
 *  @param part_type_a  particle type for which the interaction is defined
 *  @param part_type_b  particle type for which the interaction is defined
 *  @param min          @copybrief TabulatedPotential::minval
 *  @param max          @copybrief TabulatedPotential::maxval
 *  @param energy       @copybrief TabulatedPotential::energy_tab
 *  @param force        @copybrief TabulatedPotential::force_tab
 *  @retval ES_OK
 */
int tabulated_set_params(int part_type_a, int part_type_b, double min,
                         double max, std::vector<double> const &energy,
                         std::vector<double> const &force);

/** Calculate a non-bonded pair force factor by linear interpolation from a
 *  table.
 */
inline double tabulated_pair_force_factor(IA_parameters const &ia_params,
                                          double dist) {
  if (dist < ia_params.tab.cutoff()) {
    return ia_params.tab.force(dist) / dist;
  }
  return 0.0;
}

/** Calculate a non-bonded pair force by linear interpolation from a table. */
inline Utils::Vector3d tabulated_pair_force(IA_parameters const &ia_params,
                                            Utils::Vector3d const &d,
                                            double dist) {
  return d * tabulated_pair_force_factor(ia_params, dist);
}

/** Calculate a non-bonded pair energy by linear interpolation from a table. */
inline double tabulated_pair_energy(IA_parameters const &ia_params,
                                    double dist) {
  if (dist < ia_params.tab.cutoff()) {
    return ia_params.tab.energy(dist);
  }
  return 0.0;
}

#endif

#endif
