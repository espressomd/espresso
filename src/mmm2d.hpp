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
/** \file mmm2d.hpp MMM2D algorithm for long range coulomb interaction
    in 2d+h geometries.  Implementation of the MMM2D method for the
    calculation of the electrostatic interaction for two dimensionally
    periodic systems. For details on the method see MMM general. The
    MMM2D method works only with the layered or nsquared \ref
    tclcommand_cellsystem "cell system". The tuning is not automated,
    since the only tunable parameter is the cell size, which can be
    changed easily in Tcl. Moreover, only a few values make sense to
    be tested, since in general the number of cells will be between 5
    and 20.
 */
#ifndef MMM2D_H
#define MMM2D_H
#include "utils.hpp"
#ifdef ELECTROSTATICS

/** error messages, see above */
extern char const *mmm2d_errors[];

/** parameters for the MMM2D method for electrostatics. */
typedef struct {
  /** maximal error of a pairwise interaction. Used at least by the
      near formula, since this does the error control at run time */
  double maxPWerror;
  /** far formula cutoff and its square */
  double far_cut, far_cut2;
  /** if nonzero, \ref MMM2D_struct::far_cut has been calculated by \ref MMM2D_tune_far, and will be
      recalculated automatically, if important parameters, such as the number of cells, change. If this
      is zero, the far cutoff has been set explicitly by the user and will not be touched by Espresso. */
  int far_calculated;
  /// whether there is dielectric contrast
  int dielectric_contrast_on;
  /** dielectric contrasts at the bottom and top of the simulation cell */
  double delta_mid_top, delta_mid_bot, delta_mult;
} MMM2D_struct;
extern MMM2D_struct mmm2d_params;

/** set parameters for MMM2D. This assumes that the particles do NOT leave the box.
    For the near formula (nsquared cell structure), precision might be lost, while
    the far formula might have problems with overflows. 
    @param maxPWerror   the maximal error for the pairwise interactions. Both for
                        potential and force components. The potential is therefore
                        always slightly more precise
    @param far_cut      sets the cutoff for the far formula in inverse lengths.
                        If -1, the far cutoff is determined by maxPWerror.
			Manual setting is probably only good for testing
    @param delta_top    dielectric contrast at top of the simulation box
    @param delta_mid    dielectric contrast in the middle of the simulation box
*/
int MMM2D_set_params(double maxPWerror, double far_cut, double delta_top, double delta_mid);

/** the general long range force/energy calculation */
double MMM2D_add_far(int f, int e);

/** the actual long range force calculation */
inline void MMM2D_add_far_force() {
  MMM2D_add_far(1,0);
}

/** the actual long range energy calculation */
inline double MMM2D_far_energy() {
  return MMM2D_add_far(0,1);
}

/** pairwise calculated parts of MMM2D force (near neighbors) */
void add_mmm2d_coulomb_pair_force(double charge_factor,
				  double dv[3], double d2, double d, double f[3]);

/** pairwise calculated parts of MMM2D force (near neighbors) */
double mmm2d_coulomb_pair_energy(double charge_factor,
				 double dv[3], double d2, double d);

/// check that MMM2D can run with the current parameters
int MMM2D_sanity_checks();

/// initialize the MMM2D constants 
void MMM2D_init();

/** if the number of particles has changed (even per node),
    the particle buffers for the coefficients have to be resized. */
void MMM2D_on_resort_particles();

/** energy contribution from dielectric layers */
double  MMM2D_dielectric_layers_energy_contribution();

/** force contribution from dielectric layers */
void  MMM2D_dielectric_layers_force_contribution();

#endif

#endif
