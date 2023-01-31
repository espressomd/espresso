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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_WALBERLACHECKPOINT_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_WALBERLACHECKPOINT_HPP

namespace ScriptInterface::walberla {

enum class CptMode : int {
  ascii = 0,
  binary = 1,
  unit_test_runtime_error = -1,
  unit_test_ios_failure = -2
};

/** Inject code for unit tests. */
inline void unit_test_handle(int mode) {
  switch (mode) {
  case static_cast<int>(CptMode::ascii):
  case static_cast<int>(CptMode::binary):
    return;
  case static_cast<int>(CptMode::unit_test_runtime_error):
    throw std::runtime_error("unit test error");
  case static_cast<int>(CptMode::unit_test_ios_failure):
    throw std::ios_base::failure("unit test error");
  default:
    throw std::domain_error("Unknown mode " + std::to_string(mode));
  }
}

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_WALBERLACHECKPOINT_HPP
