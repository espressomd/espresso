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
#ifndef BUCKINGHAM_H
#define BUCKINGHAM_H
/** \file
 *  Routines to calculate the Buckingham potential between particle pairs.
 *
 *  Implementation in \ref buckingham.cpp.
 */

#include "config.hpp"

#ifdef BUCKINGHAM

#include "debug.hpp"
#include "nonbonded_interaction_data.hpp"

int buckingham_set_params(int part_type_a, int part_type_b, double A, double B,
                          double C, double D, double cut, double discont,
                          double shift);

/**Resultant Force due to a Buckingham potential between two particles at
 * interatomic separation r greater than or equal to discont*/
inline double buck_force_r(double A, double B, double C, double D, double r) {
  return (A * B * exp(-B * r) - 6.0 * C / pow(r, 7) - 4.0 * D / pow(r, 5));
}
/**Potential Energy due to a Buckingham potential between two particles at
 * interatomic separation r greater than or equal to discont*/
inline double buck_energy_r(double A, double B, double C, double D,
                            double shift, double r) {
  return (A * exp(-B * r) - C / pow(r, 6) - D / pow(r, 4) + shift);
}

/** Calculate Buckingham force between particle p1 and p2 and add
    it to their force. */
inline void add_buck_pair_force(Particle const *const p1,
                                Particle const *const p2,
                                IA_parameters const *const ia_params,
                                Utils::Vector3d const &d, double dist,
                                Utils::Vector3d &force) {
  if (dist < ia_params->buckingham.cut) {
    /* case: resulting force/energy greater than discontinuity and
             less than cutoff (true Buckingham region) */
    double fac;
    if (dist > ia_params->buckingham.discont) {
      fac =
          buck_force_r(ia_params->buckingham.A, ia_params->buckingham.B,
                       ia_params->buckingham.C, ia_params->buckingham.D, dist) /
          dist;
      force += fac * d;
#ifdef LJ_WARN_WHEN_CLOSE
      if (fac * dist > 1000)
        fprintf(stderr, "%d: BUCK-Warning: Pair (%d-%d) force=%f dist=%f\n",
                this_node, p1->p.identity, p2->p.identity, fac * dist, dist);
#endif
    } else {
      /* resulting force/energy in the linear region*/
      fac = -ia_params->buckingham.F2 / dist;
      force += fac * d;
    }

    ONEPART_TRACE(if (p1->p.identity == check_id)
                      fprintf(stderr,
                              "%d: OPT: BUCK   f = (%.3e,%.3e,%.3e) "
                              "with part id=%d at dist %f fac %.3e\n",
                              this_node, p1->f.f[0], p1->f.f[1], p1->f.f[2],
                              p2->p.identity, dist, fac));
    ONEPART_TRACE(if (p2->p.identity == check_id)
                      fprintf(stderr,
                              "%d: OPT: BUCK   f = (%.3e,%.3e,%.3e) "
                              "with part id=%d at dist %f fac %.3e\n",
                              this_node, p2->f.f[0], p2->f.f[1], p2->f.f[2],
                              p1->p.identity, dist, fac));

    BUCK_TRACE(fprintf(
        stderr, "%d: BUCK: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
        this_node, p1->p.identity, p2->p.identity, dist, fac * d[0], fac * d[1],
        fac * d[2]));
  }
}

/** calculate Buckingham energy between particle p1 and p2. */
inline double buck_pair_energy(Particle const *const p1,
                               Particle const *const p2,
                               IA_parameters const *const ia_params,
                               Utils::Vector3d const &d, double dist) {
  if (dist < ia_params->buckingham.cut) {
    /* case: resulting force/energy greater than discont and
             less than cutoff (true Buckingham region) */
    if (dist > ia_params->buckingham.discont) {
      return buck_energy_r(ia_params->buckingham.A, ia_params->buckingham.B,
                           ia_params->buckingham.C, ia_params->buckingham.D,
                           ia_params->buckingham.shift, dist);
    }

    /* resulting force/energy in the linear region*/
    return ia_params->buckingham.F1 + ia_params->buckingham.F2 * dist;
  }
  return 0.0;
}

#endif
#endif
