// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file mmm1d.h  MMM1D algorithm for long range coulomb interactions.
    Implementation of the MMM1D method for the calculation of the electrostatic interaction in one dimensionally
    periodic systems. For details on the method see \ref MMM_general. The MMM1D method works only with the
    nsquared \ref tcl_cellsystem "cell system", since neither the near nor far formula can be decomposed. However,
    this implementation is reasonably fast, so that one can use up to 200 charges easily in a simulation.
*/
#ifndef MMM1D_H
#define MMM1D_H

#include "utils.h"
#include "particle_data.h"

#ifdef ELECTROSTATICS

/** parameters for the MMM1D electrostatic interaction */
typedef struct {
  /** square of the switching radius */
  double far_switch_radius_2;
  /** cutoff of the bessel sum */
  int    bessel_cutoff;
  /** Wether to recalculate the Bessel cutoff automatically.
      If some parameters like the box dimensions change, the current
      Bessel cutoff may not be suitable to achieve the required accuracy
      anymore. If the user did not specify the Bessel cutoff explicitly,
      this flag is 1, and whenever necessary, the Bessel cutoff is
      recalculated.
  */
  int    bessel_calculated;
  /** required accuracy */
  double maxPWerror;
} MMM1D_struct;
extern MMM1D_struct mmm1d_params;

/// print the mmm1d parameters to the interpreters result
int printMMM1DToResult(Tcl_Interp *interp);

/// parse the mmm1d parameters
int inter_parse_mmm1d(Tcl_Interp *interp, int argc, char **argv);

/** parameters for MMM1D. Most of the parameters can also be tuned automatically. Unlike
    P3M, this tuning is redone automatically whenever parameters change, but not immediately
    if you set this parameters.
    @param switch_rad at which xy-distance the calculation switches from the far to the
                      near formula. If -1, this parameter will be tuned automatically.
    @param bessel_cutoff the cutoff for the bessel sum, aka far formula. Normally set this
                         to -1, then the cutoff is automatically determined using the error formula.
    @param maxPWerror the maximal allowed error for the potential and the forces without the
                      prefactors, i. e. for the pure lattice 1/r-sum. */
int MMM1D_set_params(double switch_rad, int bessel_cutoff, double maxPWerror);

/** tuning of the parameters which are not set by the user, e.g. the switching radius or the
    bessel_cutoff. */
int MMM1D_tune(Tcl_Interp *interp);

/** recalculate the polygamma taylor series. */
void MMM1D_recalcTables();

/// check that MMM1D can run with the current parameters
int MMM1D_sanity_checks();

/// initialize the MMM1D constants
void MMM1D_init();

///
void add_mmm1d_coulomb_pair_force(Particle *p1, Particle *p2, double d[3], double dist2,
				  double dist, double force[3]);

///
double mmm1d_coulomb_pair_energy(Particle *p1, Particle *p2, double d[3], double r2, double r);

#endif
#endif
