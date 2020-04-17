/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** \file
 *  main header-file for MDLC (magnetic dipolar layer correction).
 *
 *  Developer: Joan J. Cerda.
 *
 *  Purpose:   get the corrections for dipolar 3D algorithms
 *             when applied to a slab geometry and dipolar
 *             particles. DLC & co
 *
 *  Article:   @cite brodka04a
 *
 *  We also include a tuning function that returns the
 *  cut-off necessary to attend a certain accuracy.
 *
 *  Restrictions: the slab must be such that the z is the short
 *                direction. Otherwise we get trash.
 *
 *  Limitations:  at this moment it is restricted to work with 1 cpu
 */

#ifndef _DLC_DIPOLAR_H
#define _DLC_DIPOLAR_H

#include "config.hpp"
#include <ParticleRange.hpp>

#ifdef DP3M

/** parameters for the MDLC method */
struct DLC_struct {
  /** maximal pairwise error of the potential and force */
  double maxPWerror;

  /** Cutoff of the exponential sum. Since in all other MMM methods this is
   *  the far formula, we call it here the same, although in the ELC context
   *  it does not make much sense.
   */
  double far_cut;

  /** Size of the empty gap. Note that MDLC relies on the user to make sure
   *  that this condition is fulfilled.
   */
  double gap_size;

  /** Flag whether #far_cut was set by the user, or calculated by ESPResSo.
   *  In the latter case, the cutoff will be adapted if important parameters,
   *  such as the box dimensions, change.
   */
  int far_calculated;

  /** Up to where particles can be found */
  double h;

  template <class Archive> void serialize(Archive &ar, long int) {
    ar &maxPWerror &far_cut &gap_size &far_calculated &h;
  }
};
extern DLC_struct dlc_params;

int mdlc_set_params(double maxPWerror, double gap_size, double far_cut);
int mdlc_sanity_checks();
void add_mdlc_force_corrections(const ParticleRange &particles);
double add_mdlc_energy_corrections(const ParticleRange &particles);
/** Calculate the maximal dipole moment in the system */
void calc_mu_max();
#endif // DP3M

#endif
