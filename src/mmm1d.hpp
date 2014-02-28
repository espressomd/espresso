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
/** \file mmm1d.hpp MMM1D algorithm for long range coulomb interactions.
    Implementation of the MMM1D method for the calculation of the
    electrostatic interaction in one dimensionally periodic
    systems. For details on the method see MMM in general. The MMM1D
    method works only with the nsquared \ref tclcommand_cellsystem
    "cell system", since neither the near nor far formula can be
    decomposed. However, this implementation is reasonably fast, so
    that one can use up to 200 charges easily in a simulation.  */
#ifndef MMM1D_H
#define MMM1D_H

#include "utils.hpp"
#include "particle_data.hpp"

#ifdef ELECTROSTATICS

/** parameters for the MMM1D electrostatic interaction */
typedef struct {
  /** square of the switching radius */
  double far_switch_radius_2;
  /** required accuracy */
  double maxPWerror;
} MMM1D_struct;
extern MMM1D_struct mmm1d_params;

/** parameters for MMM1D. Most of the parameters can also be tuned automatically. Unlike
    P3M, this tuning is redone automatically whenever parameters change, but not immediately
    if you set this parameters.
    @param switch_rad at which xy-distance the calculation switches from the far to the
                      near formula. If -1, this parameter will be tuned automatically.
    @param maxPWerror the maximal allowed error for the potential and the forces without the
                      prefactors, i. e. for the pure lattice 1/r-sum. */
int MMM1D_set_params(double switch_rad, double maxPWerror);

/// check that MMM1D can run with the current parameters
int MMM1D_sanity_checks();

/// initialize the MMM1D constants
void MMM1D_init();

///
void add_mmm1d_coulomb_pair_force(double chprf, double d[3], double dist2,
				  double dist, double force[3]);

///
double mmm1d_coulomb_pair_energy(Particle *p1, Particle *p2, double d[3], double r2, double r);

/** tuning of the parameters which are not set by the user, e.g. the
    switching radius or the bessel_cutoff. Call this only on the master node.

    @param log contains information about the tuning (tried values and errors)
    @return \ref ES_OK or \ref ES_ERROR
*/
int mmm1d_tune(char **log);

#endif
#endif
