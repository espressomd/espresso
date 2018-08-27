/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef _THOLE_H
#define _THOLE_H
/** \file thole.hpp
 *  Routines to calculate the thole dampingi energy and/or force
 *  for a particle pair.
 *  \ref forces.cpp
 */

#include "config.hpp"

#ifdef THOLE

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "debug.hpp"
#include "electrostatics_magnetostatics/p3m.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

int thole_set_params(int part_type_a, int part_type_b, double scaling_coeff,
                     double q1q2);

inline void add_thole_pair_force(const Particle *const p1,
                                 const Particle *const p2,
                                 IA_parameters *ia_params, double d[3],
                                 double dist, double force[3]) {
  double thole_q1q2 = ia_params->THOLE_q1q2;
  double thole_s = ia_params->THOLE_scaling_coeff;

  if (thole_s != 0 && thole_q1q2 != 0 && dist < p3m.params.r_cut &&
      !(pair_bond_enum_exists_between(p1, p2, BONDED_IA_THERMALIZED_DIST))) {

    double dist2 = dist * dist;

    // Subtract p3m shortrange (dipole-dipole charge portion) to add damped
    // coulomb later
    p3m_add_pair_force(-thole_q1q2, d, dist2, dist, force);

    // Calc damping function (see doi.org/10.1016/0301-0104(81)85176-2)
    // S(r) = 1.0 - (1.0 + thole_s*r/2.0) * exp(-thole_s*r);
    // Calc F = - d/dr ( S(r)*q1q2/r) =
    // -(1/2)*(-2+(r^2*s^2+2*r*s+2)*exp(-s*r))*q1q2/r^2 Everything before
    // q1q2/r^2 can be used as a factor for the p3m_add_pair_force method
    double sr = thole_s * dist;
    double dS_r = 0.5 * (2.0 - (exp(-sr) * (sr * (sr + 2.0) + 2.0)));
    // Add damped p3m shortrange of dipole term
    p3m_add_pair_force(thole_q1q2 * dS_r, d, dist2, dist, force);

    ONEPART_TRACE(if (p1->p.identity == check_id)
                      fprintf(stderr,
                              "%d: OPT: THOLE   f = (%.3e,%.3e,%.3e) with part "
                              "id=%d at dist %f fac %.3e\n",
                              this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                              p2->p.identity, dist, thole_s));
    ONEPART_TRACE(if (p2->p.identity == check_id)
                      fprintf(stderr,
                              "%d: OPT: THOLE   f = (%.3e,%.3e,%.3e) with part "
                              "id=%d at dist %f fac %.3e\n",
                              this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                              p1->p.identity, dist, thole_s));
  }
}

inline double thole_pair_energy(const Particle *p1, const Particle *p2,
                                const IA_parameters *ia_params,
                                const double d[3], double dist) {
  double e_thole = 0;

  double thole_s = ia_params->THOLE_scaling_coeff;
  double thole_q1q2 = ia_params->THOLE_q1q2;

  if (thole_s != 0 && thole_q1q2 != 0 && dist < p3m.params.r_cut &&
      !(pair_bond_enum_exists_between(p1, p2, BONDED_IA_THERMALIZED_DIST))) {

    double dist2 = dist * dist;
    double chgfac = p1->p.q * p2->p.q;

    // Subtract p3m shortrange energy
    e_thole += p3m_pair_energy(-chgfac, dist);

    // Add damped p3m shortrange energy
    double sd = thole_s * dist;
    double S_r = 1.0 - (1.0 + sd / 2.0) * exp(-sd);
    e_thole += p3m_pair_energy(thole_q1q2 * S_r, dist);
  }
  return e_thole;
}

#endif
#endif
