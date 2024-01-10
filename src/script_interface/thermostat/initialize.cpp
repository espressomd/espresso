/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "initialize.hpp"

#include "thermostat.hpp"

namespace ScriptInterface {
namespace Thermostat {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<Thermostat>("Thermostat::Thermostat");
  om->register_new<Langevin>("Thermostat::Langevin");
  om->register_new<Brownian>("Thermostat::Brownian");
#ifdef NPT
  om->register_new<IsotropicNpt>("Thermostat::IsotropicNpt");
#endif
#ifdef WALBERLA
  om->register_new<LBThermostat>("Thermostat::LB");
#endif
#ifdef DPD
  om->register_new<DPDThermostat>("Thermostat::DPD");
#endif
#ifdef STOKESIAN_DYNAMICS
  om->register_new<Stokesian>("Thermostat::Stokesian");
#endif
  om->register_new<ThermalizedBond>("Thermostat::ThermalizedBond");
}

} // namespace Thermostat
} // namespace ScriptInterface
