/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef umbrella_H
#define umbrella_H

/** \file
 *  Routines to calculate the umbrella potential (harmonic interaction in only
 *  one Cartesian direction) between particle pairs. Useful for umbrella
 *  sampling simulations.
 *
 *  Implementation in \ref umbrella.cpp.
 */

#include "bonded_interactions/bonded_interaction_data.hpp"

#ifdef UMBRELLA

/** Set the parameters of an umbrella bond
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int umbrella_set_params(int bond_type, double k, int dir, double r);

/** Resultant force due to an umbrella potential */
inline double umbrella_force_r(double k, int dir, double r, double distn) {
  return -k * (distn - r);
}

/** Compute the umbrella bond length force.
 *  @param[in]  ia_params Bonded parameters for the pair interaction.
 *  @param[in]  d         %Distance between the particles.
 */
inline boost::optional<Utils::Vector3d>
umbrella_pair_force(Bonded_ia_parameters const &ia_params,
                    Utils::Vector3d const &d) {
  auto const distn = d[ia_params.p.umbrella.dir];
  auto const fac = -ia_params.p.umbrella.k * (distn - ia_params.p.umbrella.r);
  Utils::Vector3d force{};
  force[ia_params.p.umbrella.dir] = fac;

  return force;
}

/** Compute the umbrella bond length energy.
 *  @param[in]  ia_params Bonded parameters for the pair interaction.
 *  @param[in]  d         %Distance between the particles.
 */
inline boost::optional<double>
umbrella_pair_energy(Bonded_ia_parameters const &ia_params,
                     Utils::Vector3d const &d) {
  auto const distn = d[ia_params.p.umbrella.dir];
  return 0.5 * ia_params.p.umbrella.k *
         Utils::sqr(distn - ia_params.p.umbrella.r);
}

#endif /* ifdef UMBRELLA */
#endif
