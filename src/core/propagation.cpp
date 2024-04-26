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

#include "propagation.hpp"

#include "config/config.hpp"

#include "PropagationMode.hpp"

#include <cassert>
#include <string>
#include <unordered_map>

// no-op function to generate a trace for code coverage
static bool force_code_coverage(bool value) { return value; }

/**
 * Note for developers: when enabling new propagation mode combinations,
 * make sure every single line of this function has code coverage.
 */
bool is_valid_propagation_combination(int propagation) {
  using namespace PropagationMode;
  assert(propagation >= 0);
  // check allowlist
  switch (propagation) {
  // only one mode
  case NONE: // NOLINT(bugprone-branch-clone)
    return force_code_coverage(true);
  case SYSTEM_DEFAULT:
    return force_code_coverage(true);
  case TRANS_NEWTON:
    return force_code_coverage(true);
  case TRANS_BROWNIAN:
    return force_code_coverage(true);
  case TRANS_LANGEVIN:
    return force_code_coverage(true);
#ifdef NPT
  case TRANS_LANGEVIN_NPT:
    return force_code_coverage(true);
#endif
#ifdef ROTATION
  case ROT_EULER:
    return force_code_coverage(true);
  case ROT_BROWNIAN:
    return force_code_coverage(true);
  case ROT_LANGEVIN:
    return force_code_coverage(true);
#endif // ROTATION
  case TRANS_LB_MOMENTUM_EXCHANGE:
    return force_code_coverage(true);
#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
  case TRANS_LB_TRACER:
    return force_code_coverage(true);
#endif
#ifdef ROTATION
  // same mode for translation and rotation
  case TRANS_NEWTON | ROT_EULER:
    return force_code_coverage(true);
  case TRANS_LANGEVIN | ROT_LANGEVIN:
    return force_code_coverage(true);
  case TRANS_BROWNIAN | ROT_BROWNIAN:
    return force_code_coverage(true);
#ifdef STOKESIAN_DYNAMICS
  case TRANS_STOKESIAN | ROT_STOKESIAN:
    return force_code_coverage(true);
#endif
    // other mode combinations
  case TRANS_LB_MOMENTUM_EXCHANGE | ROT_EULER:
    return force_code_coverage(true);
  case TRANS_LB_MOMENTUM_EXCHANGE | ROT_LANGEVIN:
    return force_code_coverage(true);
#ifdef VIRTUAL_SITES_RELATIVE
  case TRANS_VS_RELATIVE | ROT_VS_RELATIVE:
    return force_code_coverage(true);
  case TRANS_VS_RELATIVE | ROT_VS_RELATIVE | TRANS_LB_MOMENTUM_EXCHANGE:
    return force_code_coverage(true);
  case TRANS_VS_RELATIVE | ROT_VS_RELATIVE | TRANS_LANGEVIN | ROT_LANGEVIN:
    return force_code_coverage(true);
  case TRANS_VS_RELATIVE | ROT_VS_RELATIVE | TRANS_LB_MOMENTUM_EXCHANGE |
      ROT_LANGEVIN:
    return force_code_coverage(true);
#endif // VIRTUAL_SITES_RELATIVE
#endif // ROTATION
  }
  return force_code_coverage(false);
}

std::unordered_map<std::string, int> propagation_flags_map() {
  using namespace PropagationMode;
  std::unordered_map<std::string, int> enum_values{};
  enum_values["NONE"] = NONE;
  enum_values["SYSTEM_DEFAULT"] = SYSTEM_DEFAULT;
  enum_values["TRANS_NEWTON"] = TRANS_NEWTON;
  enum_values["TRANS_LANGEVIN"] = TRANS_LANGEVIN;
  enum_values["TRANS_LANGEVIN_NPT"] = TRANS_LANGEVIN_NPT;
  enum_values["TRANS_VS_RELATIVE"] = TRANS_VS_RELATIVE;
  enum_values["TRANS_LB_MOMENTUM_EXCHANGE"] = TRANS_LB_MOMENTUM_EXCHANGE;
  enum_values["TRANS_LB_TRACER"] = TRANS_LB_TRACER;
  enum_values["TRANS_BROWNIAN"] = TRANS_BROWNIAN;
  enum_values["TRANS_STOKESIAN"] = TRANS_STOKESIAN;
  enum_values["ROT_EULER"] = ROT_EULER;
  enum_values["ROT_LANGEVIN"] = ROT_LANGEVIN;
  enum_values["ROT_VS_RELATIVE"] = ROT_VS_RELATIVE;
  enum_values["ROT_BROWNIAN"] = ROT_BROWNIAN;
  enum_values["ROT_STOKESIAN"] = ROT_STOKESIAN;
  return enum_values;
}

std::string propagation_bitmask_to_string(int propagation) {
  std::string serialized{""};
  for (auto const &[name, flag] : propagation_flags_map()) {
    if (propagation & flag) {
      serialized += "|" + name;
    }
  }
  if (serialized.empty()) {
    serialized = "NONE";
  } else {
    serialized = serialized.substr(1);
  }
  return serialized;
}
