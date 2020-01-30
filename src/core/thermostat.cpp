/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include <boost/mpi.hpp>

#include <utils/u32_to_u64.hpp>

#include "bonded_interactions/thermalized_bond.hpp"
#include "communication.hpp"
#include "dpd.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "npt.hpp"
#include "thermostat.hpp"

int thermo_switch = THERMO_OFF;
double temperature = 0.0;
bool thermo_virtual = true;

using Thermostat::GammaType;

/**
 * @brief Register a thermostat MPI callbacks
 *
 * @param thermostat        The thermostat global variable
 * @param thermostat_enum   The thermostat enum value
 */
#define REGISTER_THERMOSTAT_CALLBACKS(thermostat, thermostat_enum)             \
  void mpi_bcast_##thermostat##_rng_counter_slave(const uint64_t counter) {    \
    thermostat.rng_counter =                                                   \
        std::make_unique<Utils::Counter<uint64_t>>(counter);                   \
  }                                                                            \
                                                                               \
  REGISTER_CALLBACK(mpi_bcast_##thermostat##_rng_counter_slave)                \
                                                                               \
  void mpi_bcast_##thermostat##_rng_counter(const uint64_t counter) {          \
    mpi_call(mpi_bcast_##thermostat##_rng_counter_slave, counter);             \
  }                                                                            \
                                                                               \
  void thermostat##_rng_counter_increment() {                                  \
    if (thermo_switch & thermostat_enum)                                       \
      thermostat.rng_counter->increment();                                     \
  }                                                                            \
                                                                               \
  bool thermostat##_is_seed_required() {                                       \
    /* Seed is required if rng is not initialized */                           \
    return thermostat.rng_counter == nullptr;                                  \
  }                                                                            \
                                                                               \
  void thermostat##_set_rng_state(const uint64_t counter) {                    \
    mpi_bcast_##thermostat##_rng_counter(counter);                             \
    thermostat.rng_counter =                                                   \
        std::make_unique<Utils::Counter<uint64_t>>(counter);                   \
  }                                                                            \
                                                                               \
  uint64_t thermostat##_get_rng_state() {                                      \
    return thermostat.rng_counter->value();                                    \
  }

LangevinThermostat langevin = {};
BrownianThermostat brownian = {};
IsotropicNptThermostat npt_iso = {};

REGISTER_THERMOSTAT_CALLBACKS(langevin, THERMO_LANGEVIN)
REGISTER_THERMOSTAT_CALLBACKS(brownian, THERMO_BROWNIAN)

void thermo_init() {
  // Init thermalized bond despite of thermostat
  if (n_thermalized_bonds) {
    thermalized_bond_init();
  }
  if (thermo_switch == THERMO_OFF) {
    return;
  }
  if (thermo_switch & THERMO_LANGEVIN)
    langevin.recalc_prefactors();
#ifdef DPD
  if (thermo_switch & THERMO_DPD)
    dpd_init();
#endif
#ifdef NPT
  if (thermo_switch & THERMO_NPT_ISO) {
    if (nptiso.piston == 0.0) {
      thermo_switch = (thermo_switch ^ THERMO_NPT_ISO);
    } else {
      npt_iso.recalc_prefactors(nptiso.piston);
    }
  }
#endif
  if (thermo_switch & THERMO_BROWNIAN)
    brownian.recalc_prefactors();
}
