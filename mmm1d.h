// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef MMM1D_H
#define MMM1D_H

#include "config.h"

#ifdef ELECTROSTATICS

typedef struct {
  double far_switch_radius_2;
  int    bessel_cutoff;
  double maxPWerror;
} MMM1D_struct;
extern MMM1D_struct mmm1d_params;

int set_mmm1d_params(Tcl_Interp *interp, double switch_rad,
		     int bessel_cutoff, double maxPWerror);

void MMM1D_recalcTables();

void MMM1D_calc_forces();

void MMM1D_init();

#endif

#endif
