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
#ifndef CORE_GENERIC_HPP
#define CORE_GENERIC_HPP

/** \file
 *  Routines to calculate the energy and/or force for particle pairs from
 *  a mathematical expression.
 *
 *  Implementation in \ref nonbonded_gen.cpp.
 *  Needs feature EXPRESSION compiled in (see \ref config.hpp).
 */

#include "config.hpp"

#ifdef EXPRESSION

#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

/** Non-Bonded potentials from mathematical expression:
 *  The pair potential and pair force are given as mathematical
 *  expressions which are parsed and evaluated.
 *  @param part_type_a particle type for which the interaction is defined
 *  @param part_type_b particle type for which the interaction is defined
 *  @param max         cutoff distance
 *  @param energy      pair potential expression
 *  @param force       pair force expression
 *  @retval ES_OK
 */
int generic_set_params(int part_type_a, int part_type_b, double max,
                       std::string const &energy, std::string const &force);

/** Add a non-bonded pair force from a mathematical expression. */
inline double generic_pair_force_factor(IA_parameters const &ia_params,
                                        double dist) {
  if (dist < ia_params.gen.cutoff()) {
    return ia_params.gen.force(dist) / dist;
  }
  return 0.0;
}

/** Add a non-bonded pair force from a mathematical expression. */
inline Utils::Vector3d generic_pair_force(IA_parameters const &ia_params,
                                          Utils::Vector3d const &d,
                                          double dist) {
  return d * generic_pair_force_factor(ia_params, dist);
}

/** Add a non-bonded pair energy from a mathematical expression. */
inline double generic_pair_energy(IA_parameters const &ia_params, double dist) {
  if (dist < ia_params.gen.cutoff()) {
    return ia_params.gen.energy(dist);
  }
  return 0.0;
}

#endif

#endif
