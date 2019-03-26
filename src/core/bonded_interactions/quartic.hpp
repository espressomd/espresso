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
 *  Routines to calculate the HARMONIC Energy or/and HARMONIC force
 *  for a particle pair.
 *  \ref forces.cpp
 */

/************************************************************/

#include "bonded_interaction_data.hpp"
#include "debug.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

/** set the parameters for the quartic potential
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int quartic_set_params(int bond_type, double k0, double k1, double r,
                       double r_cut);

/** Computes the quartic bond length force.
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] force     Force.
 *  @retval 0
 */
inline int calc_quartic_pair_force(Particle const *p1, Particle const *p2,
                                   Bonded_ia_parameters const *iaparams,
                                   double const dx[3], double force[3]) {
  int i;
  double fac;
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  double dr;

  // printf("Quartic dist2 %e, dist %e\n", dist2, dist);

  if ((iaparams->p.quartic.r_cut > 0.0) && (dist > iaparams->p.quartic.r_cut))
    return 1;

  dr = dist - iaparams->p.quartic.r;

  fac = (iaparams->p.quartic.k0 * dr + iaparams->p.quartic.k1 * dr * dr * dr) /
        dist;

  for (i = 0; i < 3; i++)
    force[i] = -fac * dx[i];

  //  printf("Quartic (%d-%d), dist %e, dx %e %e %e, dr %e, f %e %e %e\n",
  //  p1->p.identity, p2->p.identity, dist, dx[0], dx[1], dx[2], dr, force[0],
  //  force[1], force[2]);

  ONEPART_TRACE(if (p1->p.identity == check_id)
                    fprintf(stderr,
                            "%d: OPT: QUARTIC f = (%.3e,%.3e,%.3e) with part "
                            "id=%d at dist %f fac %.3e\n",
                            this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                            p2->p.identity, dist2, fac));
  ONEPART_TRACE(if (p2->p.identity == check_id)
                    fprintf(stderr,
                            "%d: OPT: QUARTIC f = (%.3e,%.3e,%.3e) with part "
                            "id=%d at dist %f fac %.3e\n",
                            this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                            p1->p.identity, dist2, fac));

  return 0;
}

/** Computes the quartic bond length energy.
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] _energy   Energy.
 *  @retval 0
 */
inline int quartic_pair_energy(Particle const *p1, Particle const *p2,
                               Bonded_ia_parameters const *iaparams,
                               double const dx[3], double *_energy) {
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);

  if ((iaparams->p.quartic.r_cut > 0.0) && (dist > iaparams->p.quartic.r_cut))
    return 1;

  double dr2 = Utils::sqr(dist - iaparams->p.quartic.r);

  *_energy = 0.5 * iaparams->p.quartic.k0 * dr2 +
             0.25 * iaparams->p.quartic.k1 * Utils::sqr(dr2);
  return 0;
}

#endif
