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
#ifndef umbrella_H
#define umbrella_H

/** \file
 *  Routines to calculate the umbrella energy and/or  force
 *  for a particle pair: harmonic interaction in only one
 *  Cartesian direction. Useful for umbrella sampling simulations.
 *  \ref forces.cpp
 */

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "particle_data.hpp"

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
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  ia_params Bonded parameters for the pair interaction.
 *  @param[in]  d         %Distance between the particles.
 *  @param[out] force     Force.
 *  @retval 0
 */
inline int calc_umbrella_pair_force(Particle const *p1, Particle const *p2,
                                    Bonded_ia_parameters const *ia_params,
                                    double const d[3], double force[3]) {
  double distn;
  double fac = 0.0;
  distn = d[ia_params->p.umbrella.dir];
  fac = -ia_params->p.umbrella.k * (distn - ia_params->p.umbrella.r);
  force[ia_params->p.umbrella.dir] += fac;

  ONEPART_TRACE(if (p1->p.identity == check_id)
                    fprintf(stderr,
                            "%d: OPT: umbrella f = (%.3e,%.3e,%.3e) with part "
                            "id=%d at dist %f fac %.3e\n",
                            this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                            p2->p.identity, distn, fac));
  ONEPART_TRACE(if (p2->p.identity == check_id)
                    fprintf(stderr,
                            "%d: OPT: umbrella f = (%.3e,%.3e,%.3e) with part "
                            "id=%d at dist %f fac %.3e\n",
                            this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                            p1->p.identity, distn, fac));
  return 0;
}

/** Compute the umbrella bond length energy.
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  ia_params Bonded parameters for the pair interaction.
 *  @param[in]  d         %Distance between the particles.
 *  @param[out] _energy   Energy.
 *  @retval 0
 */
inline int umbrella_pair_energy(Particle const *p1, Particle const *p2,
                                Bonded_ia_parameters const *ia_params,
                                double const d[3], double *_energy) {
  double distn;
  distn = d[ia_params->p.umbrella.dir];
  *_energy = 0.5 * ia_params->p.umbrella.k *
             Utils::sqr(distn - ia_params->p.umbrella.r);
  return 0;
}

#endif /* ifdef UMBRELLA */
#endif
