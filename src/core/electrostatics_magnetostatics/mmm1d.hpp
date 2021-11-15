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
 *  MMM1D algorithm for long range %Coulomb interactions.
 *  Implementation of the MMM1D method for the calculation of the electrostatic
 *  interaction in one-dimensionally periodic systems. For details on the
 *  method see MMM in general. The MMM1D method works only with the nsquared,
 *  since neither the near nor far formula can be decomposed. However, this
 *  implementation is reasonably fast, so that one can use up to 200 charges
 *  easily in a simulation.
 */
#ifndef MMM1D_H
#define MMM1D_H

#include "config.hpp"

#include <utils/Vector.hpp>

#ifdef ELECTROSTATICS

/** @brief Parameters for the MMM1D electrostatic interaction */
struct MMM1DParameters {
  /** square of the switching radius */
  double far_switch_radius_2;
  /** Maximal allowed pairwise error for the potential and force.
   *  This error ignores prefactors, i.e. this is for a pure lattice 1/r-sum.
   */
  double maxPWerror;
  /** cutoff of the Bessel sum. Only used by the GPU implementation */
  int bessel_cutoff;
};
extern MMM1DParameters mmm1d_params;

/** Set parameters for MMM1D.
 *  Most of the parameters can also be tuned automatically. Unlike P3M, this
 *  tuning is redone automatically whenever parameters change, but not
 *  immediately if you set these parameters.
 *  @param switch_rad at which xy-distance the calculation switches from the far
 *      to the near formula. If -1, this parameter will be tuned automatically.
 *  @param maxPWerror @copydoc MMM1DParameters::maxPWerror
 */
void MMM1D_set_params(double switch_rad, double maxPWerror);

/// check that MMM1D can run with the current parameters
int MMM1D_sanity_checks();

/// initialize the MMM1D constants
int MMM1D_init();

void add_mmm1d_coulomb_pair_force(double chpref, Utils::Vector3d const &d,
                                  double r, Utils::Vector3d &force);

double mmm1d_coulomb_pair_energy(double q1q2, Utils::Vector3d const &d,
                                 double r2, double r);

/** Tuning of the parameters which are not set by the user. Tune either the
 *  @ref MMM1DParameters::far_switch_radius_2 "switching radius" or the
 *  @ref MMM1DParameters::bessel_cutoff "Bessel cutoff". Call this only
 *  on the head node.
 *
 *  @param verbose output information about the tuning (tried values and errors)
 *  @param timings Number of test force calculations
 *  @retval ES_OK
 *  @retval ES_ERROR
 */
int mmm1d_tune(int timings, bool verbose);

#endif
#endif
