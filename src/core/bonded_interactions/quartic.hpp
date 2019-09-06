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
#ifndef _QUARTIC_HPP
#define _QUARTIC_HPP
/** \file
 *  Routines to calculate the quartic potential between particle pairs.
 *
 *  Implementation in \ref quartic.cpp.
 */

/************************************************************/

#include "bonded_interaction_data.hpp"
#include "debug.hpp"
#include "particle_data.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

/** Set the parameters for the quartic potential
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int quartic_set_params(int bond_type, double k0, double k1, double r,
                       double r_cut);

/** Compute the quartic bond force.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] force     Force.
 *  @return whether the bond is broken
 */
inline bool calc_quartic_pair_force(Bonded_ia_parameters const *const iaparams,
                                    Utils::Vector3d const &dx,
                                    Utils::Vector3d &force) {
  auto const dist = dx.norm();

  if ((iaparams->p.quartic.r_cut > 0.0) && (dist > iaparams->p.quartic.r_cut)) {
    return true;
  }

  auto const dr = dist - iaparams->p.quartic.r;
  auto const fac =
      (iaparams->p.quartic.k0 * dr + iaparams->p.quartic.k1 * dr * dr * dr) /
      dist;
  force = -fac * dx;

  return false;
}

/** Compute the quartic bond energy.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] _energy   Energy.
 *  @return whether the bond is broken
 */
inline bool quartic_pair_energy(Bonded_ia_parameters const *const iaparams,
                                Utils::Vector3d const &dx, double *_energy) {
  auto const dist = dx.norm();

  if ((iaparams->p.quartic.r_cut > 0.0) && (dist > iaparams->p.quartic.r_cut)) {
    return true;
  }

  double dr2 = Utils::sqr(dist - iaparams->p.quartic.r);

  *_energy = 0.5 * iaparams->p.quartic.k0 * dr2 +
             0.25 * iaparams->p.quartic.k1 * Utils::sqr(dr2);
  return false;
}

#endif
