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
#ifndef _BONDED_COULOMB_P3M_SR_HPP
#define _BONDED_COULOMB_P3M_SR_HPP
/** \file bonded_coulomb.hpp
 *  Routines to calculate the BONDED_COULOMB_P3M_SR Energy or/and
 * BONDED_COULOMB_P3M_SR force
 *  for a particle pair. This is only the shortrange part of P3M and first
 *  used to subtract certain intramolecular interactions in combination with
 * thole damping
 *  \ref forces.cpp
*/

/************************************************************/

#include "config.hpp"

#ifdef P3M

#include "debug.hpp"
#include "interaction_data.hpp"
#include "p3m.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

/// set the parameters for the bonded_coulomb potential
int bonded_coulomb_p3m_sr_set_params(int bond_type, double q1q2);

/** Computes the BONDED_COULOMB_P3M_SR pair force and adds this
    force to the particle forces (see \ref interaction_data.cpp).
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  bond parameters.
    @param dx        particle distance vector
    @param force     returns force of particle 1
    @return 0.
*/
inline int calc_bonded_coulomb_p3m_sr_pair_force(Particle *p1, Particle *p2,
                                                 Bonded_ia_parameters *iaparams,
                                                 double dx[3],
                                                 double force[3]) {
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  if (dist < p3m.params.r_cut) {
    // Set to zero because p3m adds forces
    force[0] = force[1] = force[2] = 0.;

    p3m_add_pair_force(iaparams->p.bonded_coulomb_p3m_sr.q1q2, dx, dist2, dist,
                       force);

    ONEPART_TRACE(if (p1->p.identity == check_id) fprintf(
        stderr, "%d: OPT: BONDED_COULOMB_P3M_SR f = (%.3e,%.3e,%.3e) with part "
                "id=%d at dist %f\n",
        this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2], p2->p.identity, dist2));
    ONEPART_TRACE(if (p2->p.identity == check_id) fprintf(
        stderr, "%d: OPT: BONDED_COULOMB_P3M_SR f = (%.3e,%.3e,%.3e) with part "
                "id=%d at dist %f\n",
        this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2], p1->p.identity, dist2));
  }

  return 0;
}

inline int bonded_coulomb_p3m_sr_pair_energy(Particle *p1, Particle *p2,
                                             Bonded_ia_parameters *iaparams,
                                             double dx[3], double *_energy) {
  double dist2 = sqrlen(dx);
  double dist = sqrt(dist2);
  *_energy = p3m_pair_energy(iaparams->p.bonded_coulomb_p3m_sr.q1q2, dist);

  return 0;
}

#endif

#endif
