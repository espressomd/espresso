// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef MMM2D_H
#define MMM2D_H

typedef struct {
  double maxPWerror;
  double far_cut, far_cut2;
  int far_calculated;
} MMM2D_struct;
extern MMM2D_struct mmm2d_params;

/** set parameters for MMM2D */
int set_mmm2d_params(Tcl_Interp *interp, double maxPWerror, double far_cut);

/** the actual long range force calculation */
void MMM2D_add_far_force();

/** the actual long range energy calculation */
double MMM2D_far_energy();

/** pairwise calculated parts of MMM2D force (near neighbors) */
void add_mmm2d_coulomb_pair_force(Particle *p1, Particle *p2,
				  double dv[3], double d2, double d);

/** pairwise calculated parts of MMM2D force (near neighbors) */
double mmm2d_coulomb_pair_energy(Particle *p1, Particle *p2,
				 double dv[3], double d2, double d);

void MMM2D_init();

void MMM2D_on_resort_particles();

#endif
