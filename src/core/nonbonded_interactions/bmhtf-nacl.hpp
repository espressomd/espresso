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
#ifndef BMHTF_NACL_H
#define BMHTF_NACL_H
/** \file
 *  Routines to calculate the Born-Meyer-Huggins-Tosi-Fumi potential
 *  between particle pairs.
 *
 *  Implementation in \ref bmhtf-nacl.cpp.
 */

#include "nonbonded_interaction_data.hpp"
#include "particle_data.hpp"

#ifdef BMHTF_NACL

int BMHTF_set_params(int part_type_a, int part_type_b, double A, double B,
                     double C, double D, double sig, double cut);

/** Calculate BMHTF force between particle p1 and p2 */
inline void add_BMHTF_pair_force(Particle const *const p1,
                                 Particle const *const p2,
                                 IA_parameters const *const ia_params,
                                 Utils::Vector3d const &d, double dist,
                                 double dist2, Utils::Vector3d &force) {
  if (dist < ia_params->bmhtf.cut) {
    auto const pw8 = dist2 * dist2 * dist2 * dist2;
    auto const fac =
        ia_params->bmhtf.A * ia_params->bmhtf.B *
            exp(ia_params->bmhtf.B * (ia_params->bmhtf.sig - dist)) / dist -
        6 * ia_params->bmhtf.C / pw8 - 8 * ia_params->bmhtf.D / pw8 / dist2;
    force += fac * d;
  }
}

/** Calculate BMHTF potential energy between particles p1 and p2. */
inline double BMHTF_pair_energy(Particle const *const p1,
                                Particle const *const p2,
                                IA_parameters const *const ia_params,
                                Utils::Vector3d const &d, double dist,
                                double dist2) {
  if (dist < ia_params->bmhtf.cut) {
    auto const pw6 = dist2 * dist2 * dist2;
    return ia_params->bmhtf.A *
               exp(ia_params->bmhtf.B * (ia_params->bmhtf.sig - dist)) -
           ia_params->bmhtf.C / pw6 - ia_params->bmhtf.D / pw6 / dist2 +
           ia_params->bmhtf.computed_shift;
  }
  return 0.0;
}

#endif
#endif
