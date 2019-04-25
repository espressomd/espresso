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
#ifndef _BONDED_COULOMB_SR_HPP
#define _BONDED_COULOMB_SR_HPP
/** \file
 *  Routines to calculate the BONDED_COULOMB_SR Energy or/and
 * BONDED_COULOMB_SR force for a particle pair. This is only the shortrange
 * part of any coulomb interaction and first used to subtract certain
 * intramolecular interactions in combination with Thole damping \ref forces.cpp
 */

/************************************************************/

#include "config.hpp"

#ifdef ELECTROSTATICS

#include "bonded_interaction_data.hpp"
#include "debug.hpp"
#include "electrostatics_magnetostatics/coulomb_inline.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

/** set the parameters for the bonded_coulomb potential
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int bonded_coulomb_sr_set_params(int bond_type, double q1q2);

/** Computes the BONDED_COULOMB_SR pair force.
 *  @param[in]  iaparams  Interaction parameters.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] force     Force.
 *  @retval 0
 */
inline int
calc_bonded_coulomb_sr_pair_force(Bonded_ia_parameters const *iaparams,
                                  double dx[3], double force[3]) {
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  auto const forcevec =
      Coulomb::central_force(iaparams->p.bonded_coulomb_sr.q1q2, dx, dist);

  force[0] = forcevec[0];
  force[1] = forcevec[1];
  force[2] = forcevec[2];

  return 0;
}

/** Computes the BONDED_COULOMB_SR pair energy.
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  iaparams  Interaction parameters.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] _energy   Energy.
 *  @retval 0
 */
inline int bonded_coulomb_sr_pair_energy(const Particle *p1, const Particle *p2,
                                         Bonded_ia_parameters const *iaparams,
                                         double *dx, double *_energy) {
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  *_energy = Coulomb::pair_energy(p1, p2, iaparams->p.bonded_coulomb_sr.q1q2,
                                  dx, dist, dist2);
  return 0;
}

#endif

#endif
