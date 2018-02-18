/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file thermostat.cpp
    Implementation of \ref thermostat.hpp "thermostat.h"
 */
#include "thermostat.hpp"
#include "communication.hpp"
#include "lattice.hpp"
#include "npt.hpp"
#include "ghmc.hpp"
#include "lb.hpp"

#include <iostream>
#include <fstream>
#include <unistd.h>

/* thermostat switch */
int thermo_switch = THERMO_OFF;
/** Temperature */
double temperature = 0.0;

using Thermostat::GammaType;

namespace {
  /* These functions return the sentinel value for the
     langevin params, indicating that they have not been
     set yet. */
constexpr double sentinel(double) { return -1.0; }
Vector3d sentinel(Vector3d) { return {-1.0, -1.0, -1.0}; }
}

/* LANGEVIN THERMOSTAT */

/* Langevin friction coefficient gamma for translation. */
GammaType langevin_gamma = sentinel(GammaType{});
/* Friction coefficient gamma for rotation. */
GammaType langevin_gamma_rotation = sentinel(GammaType{});
GammaType langevin_pref1;
GammaType langevin_pref2;
GammaType langevin_pref2_rotation;

/* Langevin for translations */
bool langevin_trans = true;
/* Langevin for rotations */
bool langevin_rotate = true;

/* NPT ISOTROPIC THERMOSTAT */
// INSERT COMMENT
double nptiso_gamma0 = 0.0;
// INSERT COMMENT
double nptiso_gammav = 0.0;

/* GHMC THERMOSTAT */
// Number of NVE-MD steps in each GHMC cycle
int ghmc_nmd = 1;
// phi parameter for partial momentum update step in GHMC
double ghmc_phi = 0;

#ifdef MULTI_TIMESTEP
GammaType langevin_pref1_small, langevin_pref2_small;
static GammaType langevin_pref2_small_buffer;
#endif

/** buffers for the work around for the correlated random values which cool the
   system,
    and require a magical heat up whenever reentering the integrator. */
static GammaType langevin_pref2_buffer;
static GammaType langevin_pref2_rotation_buffer;

#ifdef NPT
double nptiso_pref1;
double nptiso_pref2;
double nptiso_pref3;
double nptiso_pref4;
#endif

void thermo_init_langevin() {
  langevin_pref1 = -langevin_gamma / time_step;
  langevin_pref2 = sqrt(24.0 * temperature / time_step * langevin_gamma);
  ;

#ifdef MULTI_TIMESTEP
  if (smaller_time_step > 0.) {
    langevin_pref1_small = -langevin_gamma / smaller_time_step;
#ifndef LANGEVIN_PER_PARTICLE
    langevin_pref2_small =
        sqrt(24.0 * temperature * langevin_gamma / smaller_time_step);
#endif
  } else {
    langevin_pref1_small = -langevin_gamma / time_step;
#ifndef LANGEVIN_PER_PARTICLE
    langevin_pref2_small =
        sqrt(24.0 * temperature * langevin_gamma / time_step);
#endif
  }
#endif

  /* If gamma_rotation is not set explicitly,
     use the linear one. */
  if (langevin_gamma_rotation < GammaType{}) {
    langevin_gamma_rotation = langevin_gamma;
  }

  langevin_pref2_rotation =
      sqrt(24.0 * temperature * langevin_gamma_rotation / time_step);

#ifdef PARTICLE_ANISOTROPY
#ifdef ROTATION
  THERMO_TRACE(
      fprintf(stderr, "%d: thermo_init_langevin: langevin_gamma_rotation=(%f,%f,%f), "
                      "langevin_pref2_rotation=(%f,%f,%f)",
              this_node, langevin_gamma_rotation[0], langevin_gamma_rotation[1],
              langevin_gamma_rotation[2], langevin_pref2_rotation[0],
              langevin_pref2_rotation[1], langevin_pref2_rotation[2]));
#endif
  THERMO_TRACE(fprintf(
      stderr, "%d: thermo_init_langevin: langevin_pref1=(%f,%f,%f), langevin_pref2=(%f,%f,%f)",
      this_node, langevin_pref1[0], langevin_pref1[1], langevin_pref1[2],
      langevin_pref2[0], langevin_pref2[1], langevin_pref2[2]));
#else
#ifdef ROTATION
  THERMO_TRACE(
      fprintf(stderr, "%d: thermo_init_langevin: langevin_gamma_rotation=%f, "
                      "langevin_pref2_rotation=%f",
              this_node, langevin_gamma_rotation, langevin_pref2_rotation));
#endif
  THERMO_TRACE(fprintf(
      stderr, "%d: thermo_init_langevin: langevin_pref1=%f, langevin_pref2=%f",
      this_node, langevin_pref1, langevin_pref2));
#endif
}

#ifdef NPT
void thermo_init_npt_isotropic() {
  if (nptiso.piston != 0.0) {
    nptiso_pref1 = -nptiso_gamma0 * 0.5 * time_step;
#ifdef MULTI_TIMESTEP
    if (smaller_time_step > 0.)
      nptiso_pref2 = sqrt(12.0 * temperature * nptiso_gamma0 * time_step) *
                     smaller_time_step;
    else
#endif
      nptiso_pref2 =
          sqrt(12.0 * temperature * nptiso_gamma0 * time_step) * time_step;
    nptiso_pref3 = -nptiso_gammav * (1.0 / nptiso.piston) * 0.5 * time_step;
    nptiso_pref4 = sqrt(12.0 * temperature * nptiso_gammav * time_step);
    THERMO_TRACE(fprintf(
        stderr, "%d: thermo_init_npt_isotropic: nptiso_pref1=%f, "
                "nptiso_pref2=%f, nptiso_pref3=%f, nptiso_pref4=%f \n",
        this_node, nptiso_pref1, nptiso_pref2, nptiso_pref3, nptiso_pref4));
  } else {
    thermo_switch = (thermo_switch ^ THERMO_NPT_ISO);
    THERMO_TRACE(fprintf(stderr, "%d: thermo_init_npt_isotropic: switched off "
                                 "nptiso (piston=%f; thermo_switch=%d) \n",
                         this_node, nptiso.piston, thermo_switch));
  }
}
#endif

void thermo_init() {
  if (thermo_switch == THERMO_OFF) {
    return;
  }
  if(thermo_switch & THERMO_LANGEVIN ) thermo_init_langevin();
#ifdef DPD
  if(thermo_switch & THERMO_DPD) dpd_init();
#endif
#ifdef NPT
  if (thermo_switch & THERMO_NPT_ISO)
    thermo_init_npt_isotropic();
#endif
#ifdef GHMC
  if (thermo_switch & THERMO_GHMC)
    thermo_init_ghmc();
#endif
}

void langevin_heat_up() {
  langevin_pref2_buffer = langevin_pref2;
  langevin_pref2 *= sqrt(3);

  langevin_pref2_rotation_buffer = langevin_pref2_rotation;
  langevin_pref2_rotation *= sqrt(3);

#ifdef MULTI_TIMESTEP
  langevin_pref2_small_buffer = langevin_pref2_small;
  langevin_pref2_small *= sqrt(3);
#endif
}

void thermo_heat_up() {
  if (thermo_switch & THERMO_LANGEVIN) {
    langevin_heat_up();
  }
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    dpd_heat_up();
  }
#endif
}

void langevin_cool_down() {
  langevin_pref2 = langevin_pref2_buffer;
  langevin_pref2_rotation = langevin_pref2_rotation_buffer;

#ifdef MULTI_TIMESTEP
  langevin_pref2_small = langevin_pref2_small_buffer;
#endif
}

void thermo_cool_down() {
  if (thermo_switch & THERMO_LANGEVIN) {
    langevin_cool_down();
  }
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    dpd_cool_down();
  }
#endif
}
