/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
    Implementation of \ref thermostat.hpp "thermostat.hpp"
 */
#include "thermostat.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "communication.hpp"
#include "dpd.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "npt.hpp"

#include "utils/u32_to_u64.hpp"
#include <boost/mpi.hpp>

#include <fstream>
#include <iostream>
#include <unistd.h>

/* thermostat switch */
int thermo_switch = THERMO_OFF;
/** Temperature */
double temperature = 0.0;

/** True if the thermostat should act on virtual particles. */
bool thermo_virtual = true;

using Thermostat::GammaType;

namespace {
/* These functions return the sentinel value for the
   Langevin params, indicating that they have not been
   set yet. */
constexpr double sentinel(double) { return -1.0; }
Utils::Vector3d sentinel(Utils::Vector3d) { return {-1.0, -1.0, -1.0}; }
} // namespace

/* LANGEVIN THERMOSTAT */

/* Langevin friction coefficient gamma for translation. */
GammaType langevin_gamma = sentinel(GammaType{});
/* Friction coefficient gamma for rotation. */
GammaType langevin_gamma_rotation = sentinel(GammaType{});
GammaType langevin_pref1;
GammaType langevin_pref2;
GammaType langevin_pref2_rotation;

/* NPT ISOTROPIC THERMOSTAT */
double nptiso_gamma0 = 0.0;
double nptiso_gammav = 0.0;

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

std::unique_ptr<Utils::Counter<uint64_t>> langevin_rng_counter;

void mpi_bcast_langevin_rng_counter_slave(const uint64_t counter) {
  langevin_rng_counter = std::make_unique<Utils::Counter<uint64_t>>(counter);
}

REGISTER_CALLBACK(mpi_bcast_langevin_rng_counter_slave)

void mpi_bcast_langevin_rng_counter(const uint64_t counter) {
  mpi_call(mpi_bcast_langevin_rng_counter_slave, counter);
}

void langevin_rng_counter_increment() {
  if (thermo_switch & THERMO_LANGEVIN)
    langevin_rng_counter->increment();
}

bool langevin_is_seed_required() {
  /* Seed is required if rng is not initialized */
  return langevin_rng_counter == nullptr;
}

void langevin_set_rng_state(const uint64_t counter) {
  mpi_bcast_langevin_rng_counter(counter);
  langevin_rng_counter = std::make_unique<Utils::Counter<uint64_t>>(counter);
}

uint64_t langevin_get_rng_state() { return langevin_rng_counter->value(); }

void thermo_init_langevin() {
  langevin_pref1 = -langevin_gamma;
  langevin_pref2 = sqrt(24.0 * temperature / time_step * langevin_gamma);

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
      fprintf(stderr,
              "%d: thermo_init_langevin: langevin_gamma_rotation=(%f,%f,%f), "
              "langevin_pref2_rotation=(%f,%f,%f)",
              this_node, langevin_gamma_rotation[0], langevin_gamma_rotation[1],
              langevin_gamma_rotation[2], langevin_pref2_rotation[0],
              langevin_pref2_rotation[1], langevin_pref2_rotation[2]));
#endif
  THERMO_TRACE(fprintf(stderr,
                       "%d: thermo_init_langevin: langevin_pref1=(%f,%f,%f), "
                       "langevin_pref2=(%f,%f,%f)",
                       this_node, langevin_pref1[0], langevin_pref1[1],
                       langevin_pref1[2], langevin_pref2[0], langevin_pref2[1],
                       langevin_pref2[2]));
#else
#ifdef ROTATION
  THERMO_TRACE(fprintf(stderr,
                       "%d: thermo_init_langevin: langevin_gamma_rotation=%f, "
                       "langevin_pref2_rotation=%f",
                       this_node, langevin_gamma_rotation,
                       langevin_pref2_rotation));
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

    nptiso_pref2 = sqrt(12.0 * temperature * nptiso_gamma0 * time_step);
    nptiso_pref3 = -nptiso_gammav * (1.0 / nptiso.piston) * 0.5 * time_step;
    nptiso_pref4 = sqrt(12.0 * temperature * nptiso_gammav * time_step);
    THERMO_TRACE(fprintf(stderr,
                         "%d: thermo_init_npt_isotropic: nptiso_pref1=%f, "
                         "nptiso_pref2=%f, nptiso_pref3=%f, nptiso_pref4=%f \n",
                         this_node, nptiso_pref1, nptiso_pref2, nptiso_pref3,
                         nptiso_pref4));
  } else {
    thermo_switch = (thermo_switch ^ THERMO_NPT_ISO);
    THERMO_TRACE(fprintf(stderr,
                         "%d: thermo_init_npt_isotropic: switched off "
                         "nptiso (piston=%f; thermo_switch=%d) \n",
                         this_node, nptiso.piston, thermo_switch));
  }
}
#endif

void thermo_init() {

  // Init thermalized bond despite of thermostat
  if (n_thermalized_bonds) {
    thermalized_bond_init();
  }

  if (thermo_switch == THERMO_OFF) {
    return;
  }
  if (thermo_switch & THERMO_LANGEVIN)
    thermo_init_langevin();
#ifdef DPD
  if (thermo_switch & THERMO_DPD)
    dpd_init();
#endif
#ifdef NPT
  if (thermo_switch & THERMO_NPT_ISO)
    thermo_init_npt_isotropic();
#endif
}

void langevin_heat_up() {
  langevin_pref2_buffer = langevin_pref2;
  langevin_pref2 *= sqrt(3);

  langevin_pref2_rotation_buffer = langevin_pref2_rotation;
  langevin_pref2_rotation *= sqrt(3);
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
  if (n_thermalized_bonds) {
    thermalized_bond_heat_up();
  }
}

void langevin_cool_down() {
  langevin_pref2 = langevin_pref2_buffer;
  langevin_pref2_rotation = langevin_pref2_rotation_buffer;
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
  if (n_thermalized_bonds) {
    thermalized_bond_cool_down();
  }
}
