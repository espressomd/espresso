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

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/thermalized_bond.hpp"
#include "communication.hpp"
#include "dpd.hpp"
#include "errorhandling.hpp"
#include "npt.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <boost/variant.hpp>

void Thermostat::Thermostat::recalc_prefactors(double time_step) {
  if (thermalized_bond) {
    thermalized_bond->recalc_prefactors(time_step, *(get_system().bonded_ias));
  }
  if (langevin) {
    langevin->recalc_prefactors(kT, time_step);
  }
  if (brownian) {
    brownian->recalc_prefactors(kT);
  }
#ifdef DPD
  if (dpd) {
    dpd_init(kT, time_step);
  }
#endif
#ifdef NPT
  if (npt_iso) {
    npt_iso->recalc_prefactors(kT, nptiso.piston, time_step);
  }
#endif
}

void Thermostat::Thermostat::philox_counter_increment() {
  if (thermo_switch & THERMO_LANGEVIN) {
    langevin->rng_increment();
  }
  if (thermo_switch & THERMO_BROWNIAN) {
    brownian->rng_increment();
  }
#ifdef NPT
  if (thermo_switch & THERMO_NPT_ISO) {
    npt_iso->rng_increment();
  }
#endif
#ifdef DPD
  if (thermo_switch & THERMO_DPD) {
    dpd->rng_increment();
  }
#endif
#ifdef STOKESIAN_DYNAMICS
  if (thermo_switch & THERMO_SD) {
    stokesian->rng_increment();
  }
#endif
  if (thermo_switch & THERMO_BOND) {
    thermalized_bond->rng_increment();
  }
}

void Thermostat::Thermostat::lb_coupling_deactivate() {
  if (lb) {
    if (get_system().lb.is_solver_set() and ::comm_cart.rank() == 0 and
        lb->gamma > 0.) {
      runtimeWarningMsg()
          << "Recalculating forces, so the LB coupling forces are not "
             "included in the particle force the first time step. This "
             "only matters if it happens frequently during sampling.";
    }
    lb->couple_to_md = false;
  }
}

void ThermalizedBondThermostat::recalc_prefactors(
    double time_step, BondedInteractionsMap &bonded_ias) {
  for (auto &kv : bonded_ias) {
    if (auto *bond = boost::get<ThermalizedBond>(&(*kv.second))) {
      bond->recalc_prefactors(time_step);
    }
  }
}
