/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file elc.hpp ELC algorithm for long range coulomb interactions.
    Implementation of the ELC method for the calculation of the
    electrostatic interaction in two dimensional periodic systems. For
    details on the method see MMM in general. The ELC method works
    together with any three dimensional method, which in Espresso is
    for example \ref p3m.hpp "P3M", with metallic boundary conditions.  */
#ifndef _ELC_H
#define _ELC_H

#include "utils.hpp"
#include "particle_data.hpp"

#ifdef P3M

/** parameters for the ELC method */
typedef struct {
  /** maximal pairwise error of the potential and force */
  double maxPWerror;
  /** the cutoff of the exponential sum. Since in all other MMM methods this is
      the far formula, we call it here the same, although in the ELC context it
      does not make much sense. */
  double far_cut, far_cut2;
  /** size of the empty gap. Note that ELC relies on the user to make sure that
      this condition is fulfilled. */
  double gap_size;
  /** whether the cutoff was set by the user, or calculated by Espresso. In the latter case, the
      cutoff will be adapted if important parameters, such as the box dimensions, change. */
  int far_calculated;
  /** if true, use a homogenous neutralizing background for nonneutral systems. Unlike
      the 3d case, this background adds an additional force pointing towards the system
      center, so be careful with this. */
  
  /// neutralize the box by an homogeneous background
  int neutralize;

  /// flag whether there is any dielectric contrast in the system
  int dielectric_contrast_on;

  /// dielectric constants
  double di_top, di_mid, di_bot;
  /// dielectric prefactors
  double di_mid_top, di_mid_bot, di_fac;

  /** minimal distance of two charges for which the far formula is used. For plain ELC, this equals
      gap_size, but for dielectric ELC it is only 1./3. of that. */
  double minimal_dist;
  /** layer around the dielectric contrast in which we trick around */
  double space_layer;
  /** the space that is finally left */
  double space_box;
  /** up to where particles can be found */
  double h;

} ELC_struct;
extern ELC_struct elc_params;

/** set parameters for ELC.
    @param maxPWerror the required accuracy of the potential and the force. Note that this counts for the
    plain 1/r contribution alone, without the Bjerrum length and the charge prefactor.
    @param min_dist   the gap size.
    @param far_cut    the cutoff of the exponential sum. If -1, the cutoff is automatically calculated using
    the error formulas.
    @param neutralize whether to add a neutralizing background. WARNING: This background exerts forces, which
    are dependent on the simulation box; especially the gap size enters into the value of the forces.
    @param top dielectric constant of upper part
    @param mid dielectric constant of center part
    @param bottom  dielectric constant of lower part
*/
int ELC_set_params(double maxPWerror, double min_dist, double far_cut, int neutralize,
		   double top, double mid, double bottom);

/// the force calculation 
void ELC_add_force();

/// the energy calculation
double ELC_energy();

/// check the ELC parameters
int ELC_sanity_checks();

/// initialize the ELC constants
void ELC_init();

/// resize the particle buffers
void ELC_on_resort_particles();

/// pairwise contributions from the lowest and top layers to the energy
double ELC_P3M_dielectric_layers_energy_contribution(Particle *p1, Particle *p2); 
/// pairwise contributions from the lowest and top layers to the force
void   ELC_P3M_dielectric_layers_force_contribution(Particle *p1, Particle *p2,
						    double force1[3], double force2[3]); 
/// self energies of top and bottom layers with their virtual images
double ELC_P3M_dielectric_layers_energy_self();
/// forces of particles in border layers with themselves
void ELC_P3M_self_forces();

/// assign the additional, virtual charges, used only in energy.cpp
void   ELC_p3m_charge_assign_both();
/// assign the additional, virtual charges, used only in energy.cpp
void   ELC_p3m_charge_assign_image();

/// take into account the virtual charges in the charge sums, used in energy.cpp
void   ELC_P3M_modify_p3m_sums_both();
/// take into account the virtual charges in the charge sums, used in energy.cpp
void   ELC_P3M_modify_p3m_sums_image();

/// assign the additional, virtual charges, used only in energy.cpp
void   ELC_P3M_restore_p3m_sums();

#endif

#endif
