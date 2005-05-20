// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
/** \file mmm2d.h  MMM2D algorithm for long range coulomb interaction in 2d+h geometries.
    Implementation of the MMM2D method for the calculation of the electrostatic interaction for two dimensionally
    periodic systems. For details on the method see \ref MMM_general. The MMM2D method works only with the layered
    or nsquared \ref tcl_cellsystem "cell system". The tuning is not automated, since the only tunable parameter
    is the cell size, which can be changed easily in Tcl. Moreover, only a few values make sense to be tested,
    since in general the number of cells will be between 5 and 20.
 */
#ifndef MMM2D_H
#define MMM2D_H

#ifdef ELECTROSTATICS

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
} MMM2D_struct;
extern MMM2D_struct mmm2d_params;

/// print the mmm2d parameters to the interpreters result
int printMMM2DToResult(Tcl_Interp *interp);

/// parse the mmm2d parameters
int inter_parse_mmm2d(Tcl_Interp * interp, int argc, char ** argv);

/** set parameters for MMM2D. This assumes that the particles do NOT leave the box.
    For the near formula (nsquared cell structure), precision might be lost, while
    the far formula might have problems with overflows. 
    @param maxPWerror   the maximal error for the pairwise interactions. Both for
    potential and force components. The potential is therefore always slightly
    more precise
    @param far_cut      sets the cutoff for the far formula in inverse lengths.
    if -1, the far cutoff is determined by maxPWerror. Probably only good for testing
*/
int MMM2D_set_params(double maxPWerror, double far_cut);

/** the general long range force/energy calculation */
double MMM2D_add_far(int f, int e);

/** the actual long range force calculation */
MDINLINE void MMM2D_add_far_force() {
  MMM2D_add_far(1,0);
}

/** the actual long range energy calculation */
MDINLINE double MMM2D_far_energy() {
  return MMM2D_add_far(0,1);
}

/** pairwise calculated parts of MMM2D force (near neighbors) */
void add_mmm2d_coulomb_pair_force(Particle *p1, Particle *p2,
				  double dv[3], double d2, double d, double f[3]);

/** pairwise calculated parts of MMM2D force (near neighbors) */
double mmm2d_coulomb_pair_energy(Particle *p1, Particle *p2,
				 double dv[3], double d2, double d);

/// check that MMM2D can run with the current parameters
int MMM2D_sanity_checks();

/// initialize the MMM2D constants 
void MMM2D_init();

/** if the number of particles has changed (even per node),
    the particle buffers for the coefficients have to be resized. */
void MMM2D_on_resort_particles();

#endif

#endif
