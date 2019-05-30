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
#ifndef _FENE_HPP
#define _FENE_HPP
/** \file
 *  Routines to calculate the FENE Energy or/and FENE force
 *  for a particle pair.
 *  \ref forces.cpp
 */

#include "bonded_interaction_data.hpp"
#include "particle_data.hpp"

/************************************************************/

/** set the parameters for the fene potential
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int fene_set_params(int bond_type, double k, double drmax, double r0);

/** Computes the FENE bond length force.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] force     Force.
 *  @retval 1 if the bond is broken
 *  @retval 0 otherwise
 */
inline int calc_fene_pair_force(Bonded_ia_parameters const *iaparams,
                                Utils::Vector3d const &dx, double *force) {
  auto const len = dx.norm();

  const double dr = len - iaparams->p.fene.r0;

  if (dr >= iaparams->p.fene.drmax)
    return 1;

  double fac =
      -iaparams->p.fene.k * dr / ((1.0 - dr * dr * iaparams->p.fene.drmax2i));
  if (len > ROUND_ERROR_PREC) {
    fac /= len;
  } else {
    fac = 0.0;
  }

  for (int i = 0; i < 3; i++)
    force[i] = fac * dx[i];

  return 0;
}

/** Computes the FENE bond length force.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] _energy   Energy.
 *  @retval 1 if the bond is broken
 *  @retval 0 otherwise
 */
inline int fene_pair_energy(Bonded_ia_parameters const *iaparams,
                            Utils::Vector3d const &dx, double *_energy) {
  /* compute bond stretching (r-r0) */
  double const dr = dx.norm() - iaparams->p.fene.r0;

  /* check bond stretching */
  if (dr >= iaparams->p.fene.drmax) {
    return 1;
  }

  double energy = -0.5 * iaparams->p.fene.k * iaparams->p.fene.drmax2;
  energy *= log((1.0 - dr * dr * iaparams->p.fene.drmax2i));
  *_energy = energy;
  return 0;
}

#endif
