/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef DEBYE_HUECKEL_H
#define DEBYE_HUECKEL_H
/** \file debye_hueckel.hpp
 *  Routines to calculate the Debye_Hueckel  Energy or/and Debye_Hueckel force
 *  for a particle pair.
 */
#include "config.hpp"

#ifdef ELECTROSTATICS

#include "interaction_data.hpp"
#include "debug.hpp"

/** Structure to hold Debye-Hueckel Parameters. */
typedef struct {
  /** Cutoff for Debey-Hueckel interaction. */
  double r_cut;
  /** Debye kappa (inverse Debye length) . */
  double kappa;
#ifdef COULOMB_DEBYE_HUECKEL
  double eps_int, eps_ext;
  /** Transition distances **/
  double r0, r1;
  double alpha;
#endif
} Debye_hueckel_params;

/** Structure containing the Debye-Hueckel parameters. */
extern Debye_hueckel_params dh_params;

/** \name Functions */
/************************************************************/
/*@{*/

int dh_set_params(double kappa, double r_cut);
int dh_set_params_cdh(double kappa, double r_cut, double eps_int,
                      double eps_ext, double r0, double r1, double alpha);

/** Computes the Debye_Hueckel pair force and adds this
    force to the particle forces.
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param d         Vector pointing from p1 to p2.
    @param dist      Distance between p1 and p2.
    @param force     returns the force on particle 1.
*/
#ifdef COULOMB_DEBYE_HUECKEL
inline void add_dh_coulomb_pair_force(Particle *p1, Particle *p2, double d[3],
                                      double dist, double force[3]) {
  double fac;

  const double q1 = p1->p.q;
  const double q2 = p2->p.q;

  if ((q1 == 0.0) || (q2 == 0.0))
    return;

  if (dist >= dh_params.r_cut)
    return;

  if (dist > dh_params.r1) {
    const double kappa_dist = dh_params.kappa * dist;
    fac = coulomb.prefactor * q1 * q2 * exp(-kappa_dist) /
          (dh_params.eps_ext * dist * dist * dist) * (1.0 + kappa_dist);
  } else if (dist > dh_params.r0)
    fac = coulomb.prefactor * q1 * q2 *
          exp(dh_params.alpha * (dh_params.r0 - dist)) /
          (dh_params.eps_int * dist * dist * dist) *
          (1. - dh_params.alpha * dist);
  else
    fac = coulomb.prefactor * q1 * q2 * p2->p.q /
          (dh_params.eps_int * dist * dist * dist);

  force[0] += fac * d[0];
  force[1] += fac * d[1];
  force[2] += fac * d[2];
}
#else
inline void add_dh_coulomb_pair_force(Particle *p1, Particle *p2, double d[3],
                                      double dist, double force[3]) {
  if (dist < dh_params.r_cut) {
    double fac;
    if (dh_params.kappa > 0.0) {
      /* debye hueckel case: */
      double kappa_dist = dh_params.kappa * dist;
      fac = coulomb.prefactor * p1->p.q * p2->p.q *
            (exp(-kappa_dist) / (dist * dist * dist)) * (1.0 + kappa_dist);
    } else {
      /* pure coulomb case: */
      fac = coulomb.prefactor * p1->p.q * p2->p.q / (dist * dist * dist);
    }
    for (int j = 0; j < 3; j++)
      force[j] += fac * d[j];

    ONEPART_TRACE(if (p1->p.identity == check_id)
                      fprintf(stderr, "%d: OPT: DH   f = (%.3e,%.3e,%.3e) with "
                                      "part id=%d at dist %f fac %.3e\n",
                              this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                              p2->p.identity, dist, fac));
    ONEPART_TRACE(if (p2->p.identity == check_id)
                      fprintf(stderr, "%d: OPT: DH   f = (%.3e,%.3e,%.3e) with "
                                      "part id=%d at dist %f fac %.3e\n",
                              this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                              p1->p.identity, dist, fac));
  }
}
#endif

#ifdef COULOMB_DEBYE_HUECKEL
inline double dh_coulomb_pair_energy(Particle *p1, Particle *p2, double dist) {
  const double q1 = p1->p.q;
  const double q2 = p2->p.q;

  if ((q1 == 0.0) || (q2 == 0.0))
    return 0.0;

  const double fac = q1 * q2 * coulomb.prefactor;

  if ((dist < dh_params.r_cut)) {
    if (dist < dh_params.r0) {
      return fac / (dh_params.eps_int * dist);
    }
    if (dist < dh_params.r1) {
      return fac * exp(dh_params.alpha * (dh_params.r0 - dist)) /
             (dh_params.eps_int * dist);
    }
    return fac * exp(-dh_params.kappa * dist) / (dh_params.eps_ext * dist);
  }
  return 0.0;
}
#else
inline double dh_coulomb_pair_energy(Particle *p1, Particle *p2, double dist) {
  if (dist < dh_params.r_cut) {
    if (dh_params.kappa > 0.0)
      return coulomb.prefactor * p1->p.q * p2->p.q *
             exp(-dh_params.kappa * dist) / dist;
    else
      return coulomb.prefactor * p1->p.q * p2->p.q / dist;
  }
  return 0.0;
}
#endif
/*@}*/
#endif

#endif
