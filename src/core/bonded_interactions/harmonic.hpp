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
#ifndef _HARMONIC_HPP
#define _HARMONIC_HPP
/** \file
 *  Routines to calculate the HARMONIC Energy or/and HARMONIC force
 *  for a particle pair.
 *  \ref forces.cpp
 */

/************************************************************/

#include "bonded_interaction_data.hpp"
#include "particle_data.hpp"

#include <utils/math/sqr.hpp>

/** set the parameters for the harmonic potential
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int harmonic_set_params(int bond_type, double k, double r, double r_cut);

/** Computes the harmonic bond length force.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] force     Force.
 *  @retval 1 if the bond is broken
 *  @retval 0 otherwise
 */
inline int calc_harmonic_pair_force(Bonded_ia_parameters const *iaparams,
                                    Utils::Vector3d const &dx, double *force) {
  auto const dist = dx.norm();

  if ((iaparams->p.harmonic.r_cut > 0.0) && (dist > iaparams->p.harmonic.r_cut))
    return 1;

  auto const dr = dist - iaparams->p.harmonic.r;
  auto fac = -iaparams->p.harmonic.k * dr;
  if (dist > ROUND_ERROR_PREC) { /* Regular case */
    fac /= dist;
  } else {
    fac = 0;
  }

  for (int i = 0; i < 3; i++)
    force[i] = fac * dx[i];

  return 0;
}

/** Computes the harmonic bond length energy.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] _energy   Energy.
 *  @retval 1 if the bond is broken
 *  @retval 0 otherwise
 */
inline int harmonic_pair_energy(Bonded_ia_parameters const *iaparams,
                                Utils::Vector3d const &dx, double *_energy) {
  auto const dist = dx.norm();

  if ((iaparams->p.harmonic.r_cut > 0.0) && (dist > iaparams->p.harmonic.r_cut))
    return 1;

  *_energy =
      0.5 * iaparams->p.harmonic.k * Utils::sqr(dist - iaparams->p.harmonic.r);
  return 0;
}

#endif
