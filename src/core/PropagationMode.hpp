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

#pragma once

namespace PropagationMode {
/** @brief Flags to create bitmasks for propagation modes. */
enum PropagationMode : int {
  NONE = 0,
  SYSTEM_DEFAULT = 1 << 0,
  TRANS_NEWTON = 1 << 1,
  TRANS_LANGEVIN = 1 << 2,
  TRANS_LANGEVIN_NPT = 1 << 3,
  TRANS_VS_RELATIVE = 1 << 4,
  TRANS_LB_MOMENTUM_EXCHANGE = 1 << 5,
  TRANS_LB_TRACER = 1 << 6,
  TRANS_BROWNIAN = 1 << 7,
  TRANS_STOKESIAN = 1 << 8,
  ROT_EULER = 1 << 10,
  ROT_LANGEVIN = 1 << 11,
  ROT_VS_RELATIVE = 1 << 12,
  ROT_BROWNIAN = 1 << 13,
  ROT_STOKESIAN = 1 << 14,
};
} // namespace PropagationMode

/** @brief Integrator identifier. */
enum IntegratorSwitch : int {
  INTEG_METHOD_NPT_ISO = 0,
  INTEG_METHOD_NVT = 1,
  INTEG_METHOD_STEEPEST_DESCENT = 2,
  INTEG_METHOD_BD = 3,
  INTEG_METHOD_SD = 4,
};

/** @brief Thermostat flags. */
enum ThermostatFlags : int {
  THERMO_OFF = 0,
  THERMO_LANGEVIN = 1 << 0,
  THERMO_BROWNIAN = 1 << 1,
  THERMO_NPT_ISO = 1 << 2,
  THERMO_LB = 1 << 3,
  THERMO_SD = 1 << 4,
  THERMO_DPD = 1 << 5,
  THERMO_BOND = 1 << 6,
};
