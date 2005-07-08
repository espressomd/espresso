// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
/** \file elc.h
    ELC algorithm for long range coulomb interactions.
    Implementation of the ELC method for the calculation of the electrostatic interaction in two dimensional periodic
    systems. For details on the method see \ref MMM_general. The ELC method works together with any three dimensional
    method, which in Espresso is just \ref tcl_p3m "P3M", with metallic boundary conditions.
*/
#ifndef ELC_H
#define ELC_H

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
  double minimal_dist;
  /** whether the cutoff was set by the user, or calculated by Espresso. In the latter case, the
      cutoff will be adapted if important parameters, such as the box dimensions, change. */
  int far_calculated;
  /** if true, use a homogenous neutralizing background for nonneutral systems. Unlike
      the 3d case, this background adds an additional force pointing towards the system
      center, so be careful with this. */
  int neutralize;
} ELC_struct;
extern ELC_struct elc_params;

/// print the elc parameters to the interpreters result
int printELCToResult(Tcl_Interp *interp);

/// parse the elc parameters
int inter_parse_elc_params(Tcl_Interp * interp, int argc, char ** argv);

/** set parameters for ELC.
    @param maxPWerror the required accuracy of the potential and the force. Note that this counts for the
    plain 1/r contribution alone, without the Bjerrum length and the charge prefactor.
    @param min_dist   the gap size.
    @param far_cut    the cutoff of the exponential sum. If -1, the cutoff is automatically calculated using
    the error formulas.
    @param neutralize whether to add a neutralizing background. WARNING: This background exerts forces, which
    are dependent on the simulation box; especially the gap size enters into the value of the forces.
*/
int ELC_set_params(double maxPWerror, double min_dist, double far_cut, int neutralize);

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

#endif
