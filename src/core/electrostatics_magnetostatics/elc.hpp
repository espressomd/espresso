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
/** \file
   ELC algorithm for long range
   Coulomb interactions. Implementation of the ELC method for the calculation of
   the electrostatic interaction in two dimensional periodic systems. For
    details on the method see MMM in general. The ELC method works
    together with any three dimensional method, which in Espresso is
    for example \ref p3m.hpp "P3M", with metallic boundary conditions.  */
#ifndef _ELC_H
#define _ELC_H

#include "particle_data.hpp"

#ifdef P3M

/** @brief Parameters for the ELC method */
typedef struct {
  /** @copybrief MMM2D_struct::maxPWerror */
  double maxPWerror;
  /** Cutoff of the exponential sum. Since in all other MMM methods this is
   *  the far formula, we call it here the same, although in the ELC context
   *  it does not make much sense.
   */
  double far_cut;
  /** Squared value of #far_cut. */
  double far_cut2;
  /** Size of the empty gap. Note that ELC relies on the user to make sure
   *  that this condition is fulfilled.
   */
  double gap_size;
  /** @copybrief MMM2D_struct::far_calculated */
  int far_calculated;
  /** Flag whether the box is neutralized by an homogeneous background.
   *  If true, use a homogeneous neutralizing background for nonneutral
   *  systems. Unlike the 3D case, this background adds an additional
   *  force pointing towards the system center, so be careful with this.
   */
  int neutralize;

  /// @copybrief MMM2D_struct::dielectric_contrast_on
  int dielectric_contrast_on;

  /// @copybrief MMM2D_struct::delta_mid_top
  double delta_mid_top;
  /// @copybrief MMM2D_struct::delta_mid_bot
  double delta_mid_bot;

  /// @copybrief MMM2D_struct::const_pot_on
  int const_pot;
  /// @copybrief MMM2D_struct::pot_diff
  double pot_diff;

  /** Minimal distance of two charges for which the far formula is used.
   *  For plain ELC, this equals #gap_size, but for dielectric ELC it is
   *  only 1/3 of that.
   */
  double minimal_dist;
  /** Layer around the dielectric contrast in which we trick around. */
  double space_layer;
  /** The space that is finally left. */
  double space_box;
  /** Up to where particles can be found. */
  double h;

} ELC_struct;
extern ELC_struct elc_params;

/** Set parameters for ELC.
 *  @param maxPWerror    @copybrief ELC_struct::maxPWerror
 *                       Note that this counts for the plain 1/r contribution
 *                       alone, without the prefactor and the charge prefactor.
 *  @param min_dist      @copybrief ELC_struct::minimal_dist
 *  @param far_cut       @copybrief ELC_struct::far_cut
 *                       If -1, the cutoff is automatically calculated using
 *                       the error formulas.
 *  @param neutralize    whether to add a neutralizing background.
 *                       WARNING: This background exerts forces, which are
 *                       dependent on the simulation box; especially the gap
 *                       size enters into the value of the forces.
 *  @param delta_mid_top @copybrief ELC_struct::delta_mid_top
 *  @param delta_mid_bot @copybrief ELC_struct::delta_mid_bot
 *  @param const_pot     @copybrief ELC_struct::const_pot
 *  @param pot_diff      @copybrief ELC_struct::pot_diff
 *  @retval ES_OK
 */
int ELC_set_params(double maxPWerror, double min_dist, double far_cut,
                   int neutralize, double delta_mid_top, double delta_mid_bot,
                   int const_pot, double pot_diff);

/// the force calculation
void ELC_add_force();

/// the energy calculation
double ELC_energy();

/// check the ELC parameters
/// @retval ES_OK
/// @retval ES_ERROR
int ELC_sanity_checks();

/// initialize the ELC constants
void ELC_init();

/// resize the particle buffers
void ELC_on_resort_particles();

/// pairwise contributions from the lowest and top layers to the energy
double ELC_P3M_dielectric_layers_energy_contribution(const Particle *p1,
                                                     const Particle *p2);
/// pairwise contributions from the lowest and top layers to the force
void ELC_P3M_dielectric_layers_force_contribution(const Particle *p1,
                                                  const Particle *p2,
                                                  double *force1,
                                                  double *force2);
/// self energies of top and bottom layers with their virtual images
double ELC_P3M_dielectric_layers_energy_self();
/// forces of particles in border layers with themselves
void ELC_P3M_self_forces();

/// assign the additional, virtual charges, used only in energy.cpp
void ELC_p3m_charge_assign_both();
/// assign the additional, virtual charges, used only in energy.cpp
void ELC_p3m_charge_assign_image();

/// take into account the virtual charges in the charge sums, used in energy.cpp
void ELC_P3M_modify_p3m_sums_both();
/// take into account the virtual charges in the charge sums, used in energy.cpp
void ELC_P3M_modify_p3m_sums_image();

/// assign the additional, virtual charges, used only in energy.cpp
void ELC_P3M_restore_p3m_sums();

#endif

#endif
