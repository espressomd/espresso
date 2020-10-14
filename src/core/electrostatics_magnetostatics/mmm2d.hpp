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
    MMM2D algorithm for long range Coulomb interaction
    in 2d+h geometries.  Implementation of the MMM2D method for the calculation
    of the electrostatic interaction for two-dimensionally periodic systems.
    For details on the method see MMM general. The MMM2D method works only with
    the layered or nsquared  "cell system".  The tuning is not automated, since
    the only tunable parameter is the cell size, which can be changed easily in
    the script interface. Moreover, only a few values make sense to be tested,
    since in general the number of cells will be between 5 and 20.
 */
#ifndef MMM2D_H
#define MMM2D_H

#include "config.hpp"

#include <ParticleRange.hpp>
#include <utils/Vector.hpp>

#ifdef ELECTROSTATICS

/** MMM2D error messages */
extern char const *mmm2d_errors[];

/** @brief Parameters for the MMM2D method for electrostatics. */
typedef struct {
  /** Maximal allowed pairwise error for the potential and force.
   *  Used at least by the near formula, since this does the error control at
   *  runtime.
   */
  double maxPWerror;
  /** far formula cutoff in inverse lengths. */
  double far_cut;
  /** squared value of #far_cut. */
  double far_cut2;
  /** Flag whether #far_cut was set by the user, or calculated by Espresso.
   *  In the latter case, the cutoff will be adapted if important parameters,
   *  such as the box dimensions, change.
   */
  bool far_calculated;
  /// flag whether there is any dielectric contrast in the system.
  bool dielectric_contrast_on;
  /// @brief Flag whether a const. potential is applied.
  bool const_pot;
  /// @brief Const. potential.
  double pot_diff;
  /// dielectric contrast in the upper part of the simulation cell.
  double delta_mid_top;
  /// dielectric contrast in the lower part of the simulation cell.
  double delta_mid_bot;
  /// #delta_mid_top times #delta_mid_bot
  double delta_mult;
} MMM2D_struct;
extern MMM2D_struct mmm2d_params;

/** Set parameters for MMM2D.
 *  This assumes that the particles do NOT leave the box. For the near formula
 *  (nsquared cell structure), precision might be lost, while the far formula
 *  might have problems with overflows.
 *  @param maxPWerror    @copydoc MMM2D_struct::maxPWerror
 *  @param far_cut       sets the @copybrief MMM2D_struct::far_cut cutoff.
 *                       If -1, the far cutoff is determined by @p maxPWerror.
 *                       Manual setting is probably only good for testing.
 *  @param delta_top     @copybrief MMM2D_struct::delta_mid_top
 *  @param delta_bot     @copybrief MMM2D_struct::delta_mid_bot
 *  @param const_pot  @copybrief MMM2D_struct::const_pot
 *  @param pot_diff      @copybrief MMM2D_struct::pot_diff
 */
int MMM2D_set_params(double maxPWerror, double far_cut, double delta_top,
                     double delta_bot, bool const_pot, double pot_diff);

/** the general long range force/energy calculation */
double MMM2D_add_far(bool calc_forces, bool calc_energies,
                     const ParticleRange &particles);

/** the actual long range force calculation */
inline void MMM2D_add_far_force(const ParticleRange &particles) {
  MMM2D_add_far(true, false, particles);
}

/** the actual long range energy calculation */
inline double MMM2D_far_energy(const ParticleRange &particles) {
  return MMM2D_add_far(false, true, particles);
}

/** pairwise calculated parts of MMM2D force (near neighbors) */
void add_mmm2d_coulomb_pair_force(double pref, Utils::Vector3d const &d,
                                  double dl, Utils::Vector3d &force);

/** pairwise calculated parts of MMM2D force (near neighbors) */
double mmm2d_coulomb_pair_energy(double charge_factor,
                                 Utils::Vector3d const &dv, double d);

/// check that MMM2D can run with the current parameters
int MMM2D_sanity_checks();

/// initialize the MMM2D constants
void MMM2D_init();

/** energy contribution from dielectric layers */
double MMM2D_dielectric_layers_energy_contribution();

/** force contribution from dielectric layers */
void MMM2D_dielectric_layers_force_contribution();

#endif

#endif
