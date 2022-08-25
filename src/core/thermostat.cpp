/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** \file
 *  Implementation of \ref thermostat.hpp.
 */

#include "config/config.hpp"

#include "bonded_interactions/thermalized_bond.hpp"
#include "bonded_interactions/thermalized_bond_utils.hpp"
#include "communication.hpp"
#include "dpd.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "integrate.hpp"
#include "npt.hpp"
#include "thermostat.hpp"

#include <boost/mpi.hpp>

#include <algorithm>
#include <cstdint>

int thermo_switch = THERMO_OFF;
double temperature = 0.0;
bool thermo_virtual = true;

using Thermostat::GammaType;

/**
 * @brief Create MPI callbacks of thermostat objects
 *
 * @param thermostat        The thermostat object name
 */
#define REGISTER_THERMOSTAT_CALLBACKS(thermostat)                              \
  void mpi_##thermostat##_set_rng_seed(uint32_t const seed) {                  \
    (thermostat).rng_initialize(seed);                                         \
  }                                                                            \
                                                                               \
  REGISTER_CALLBACK(mpi_##thermostat##_set_rng_seed)                           \
                                                                               \
  void thermostat##_set_rng_seed(uint32_t const seed) {                        \
    mpi_call_all(mpi_##thermostat##_set_rng_seed, seed);                       \
  }                                                                            \
                                                                               \
  void mpi_##thermostat##_set_rng_counter(uint64_t const value) {              \
    (thermostat).set_rng_counter(value);                                       \
  }                                                                            \
                                                                               \
  REGISTER_CALLBACK(mpi_##thermostat##_set_rng_counter)                        \
                                                                               \
  void thermostat##_set_rng_counter(uint64_t const value) {                    \
    mpi_call_all(mpi_##thermostat##_set_rng_counter, value);                   \
  }

LangevinThermostat langevin = {};
BrownianThermostat brownian = {};
#ifdef NPT
IsotropicNptThermostat npt_iso = {};
#endif
ThermalizedBondThermostat thermalized_bond = {};
#ifdef DPD
DPDThermostat dpd = {};
#endif
#ifdef STOKESIAN_DYNAMICS
StokesianThermostat stokesian = {};
#endif

REGISTER_THERMOSTAT_CALLBACKS(langevin)
REGISTER_THERMOSTAT_CALLBACKS(brownian)
#ifdef NPT
REGISTER_THERMOSTAT_CALLBACKS(npt_iso)
#endif
REGISTER_THERMOSTAT_CALLBACKS(thermalized_bond)
#ifdef DPD
REGISTER_THERMOSTAT_CALLBACKS(dpd)
#endif
#ifdef STOKESIAN_DYNAMICS
REGISTER_THERMOSTAT_CALLBACKS(stokesian)
#endif

void thermo_init(double time_step) {
  // initialize thermalized bond regardless of the current thermostat
  if (n_thermalized_bonds) {
    thermalized_bond_init(time_step);
  }
  if (thermo_switch == THERMO_OFF) {
    return;
  }
  if (thermo_switch & THERMO_LANGEVIN)
    langevin.recalc_prefactors(temperature, time_step);
#ifdef DPD
  if (thermo_switch & THERMO_DPD)
    dpd_init(temperature, time_step);
#endif
#ifdef NPT
  if (thermo_switch & THERMO_NPT_ISO) {
    npt_iso.recalc_prefactors(temperature, nptiso.piston, time_step);
  }
#endif
  if (thermo_switch & THERMO_BROWNIAN)
    brownian.recalc_prefactors(temperature);
}

void philox_counter_increment() {
  if (thermo_switch & THERMO_LANGEVIN) {
    langevin.rng_increment();
  }
  if (thermo_switch & THERMO_BROWNIAN) {
    brownian.rng_increment();
  }
#ifdef NPT
  if (thermo_switch & THERMO_NPT_ISO) {
    npt_iso.rng_increment();
  }
#endif
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    dpd.rng_increment();
  }
#endif
#ifdef STOKESIAN_DYNAMICS
  if (thermo_switch & THERMO_SD) {
    stokesian.rng_increment();
  }
#endif
  if (n_thermalized_bonds) {
    thermalized_bond.rng_increment();
  }
}

void mpi_set_brownian_gamma_local(GammaType const &gamma) {
  brownian.gamma = gamma;
}

void mpi_set_brownian_gamma_rot_local(GammaType const &gamma) {
  brownian.gamma_rotation = gamma;
}

void mpi_set_langevin_gamma_local(GammaType const &gamma) {
  langevin.gamma = gamma;
  on_thermostat_param_change();
}

void mpi_set_langevin_gamma_rot_local(GammaType const &gamma) {
  langevin.gamma_rotation = gamma;
  on_thermostat_param_change();
}

REGISTER_CALLBACK(mpi_set_brownian_gamma_local)
REGISTER_CALLBACK(mpi_set_brownian_gamma_rot_local)
REGISTER_CALLBACK(mpi_set_langevin_gamma_local)
REGISTER_CALLBACK(mpi_set_langevin_gamma_rot_local)

void mpi_set_brownian_gamma(GammaType const &gamma) {
  mpi_call_all(mpi_set_brownian_gamma_local, gamma);
}

void mpi_set_brownian_gamma_rot(GammaType const &gamma) {
  mpi_call_all(mpi_set_brownian_gamma_rot_local, gamma);
}

void mpi_set_langevin_gamma(GammaType const &gamma) {
  mpi_call_all(mpi_set_langevin_gamma_local, gamma);
}
void mpi_set_langevin_gamma_rot(GammaType const &gamma) {
  mpi_call_all(mpi_set_langevin_gamma_rot_local, gamma);
}

void mpi_set_thermo_virtual_local(bool thermo_virtual) {
  ::thermo_virtual = thermo_virtual;
}

REGISTER_CALLBACK(mpi_set_thermo_virtual_local)

void mpi_set_thermo_virtual(bool thermo_virtual) {
  mpi_call_all(mpi_set_thermo_virtual_local, thermo_virtual);
}

void mpi_set_temperature_local(double temperature) {
  ::temperature = temperature;
  try {
    on_temperature_change();
  } catch (std::exception const &err) {
    runtimeErrorMsg() << err.what();
  }
  on_thermostat_param_change();
}

REGISTER_CALLBACK(mpi_set_temperature_local)

void mpi_set_temperature(double temperature) {
  mpi_call_all(mpi_set_temperature_local, temperature);
}

void mpi_set_thermo_switch_local(int thermo_switch) {
  ::thermo_switch = thermo_switch;
}

REGISTER_CALLBACK(mpi_set_thermo_switch_local)

void mpi_set_thermo_switch(int thermo_switch) {
  mpi_call_all(mpi_set_thermo_switch_local, thermo_switch);
}

#ifdef NPT
void mpi_set_nptiso_gammas_local(double gamma0, double gammav) {
  npt_iso.gamma0 = gamma0;
  npt_iso.gammav = gammav;
  on_thermostat_param_change();
}

REGISTER_CALLBACK(mpi_set_nptiso_gammas_local)

void mpi_set_nptiso_gammas(double gamma0, double gammav) {
  mpi_call_all(mpi_set_nptiso_gammas_local, gamma0, gammav);
}
#endif
