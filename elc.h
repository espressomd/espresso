// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef ELC_H
#define ELC_H

typedef struct {
  double maxPWerror;
  double far_cut, far_cut2;
  double minimal_dist;
  int far_calculated;
} ELC_struct;
extern ELC_struct elc_params;

/// print the elc parameters to the interpreters result
int printELCToResult(Tcl_Interp *interp);

/// parse the elc parameters
int inter_parse_elc_params(Tcl_Interp * interp, int argc, char ** argv);

/** set parameters for ELC */
int ELC_set_params(double maxPWerror, double min_dist, double far_cut);

void ELC_add_force();

double ELC_energy();

int ELC_sanity_checks();

void ELC_init();

void ELC_on_resort_particles();

#endif
