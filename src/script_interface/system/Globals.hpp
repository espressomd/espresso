/*
 * Copyright (C) 2013-2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_SYSTEM_GLOBALS_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_SYSTEM_GLOBALS_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <cassert>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace System {

class Globals : public AutoParameters<Globals> {
public:
  Globals();

private:
  void do_construct(VariantMap const &params) override {
    /* When reloading the system state from a checkpoint file,
     * the order of global variables instantiation matters.
     * The @c box_l must be set before any other global variable.
     * All these globals re-initialize the cell system, and we
     * cannot use the default-constructed @c box_geo when e.g.
     * long-range interactions exist in the system, otherwise
     * runtime errors about the local geometry being smaller
     * than the interaction range would be raised.
     */
    if (not params.empty()) {
      auto const keys =
          std::vector<std::string>{"box_l", "periodicity", "min_global_cut"};
      assert(params.size() == keys.size());
      for (auto const &key : keys) {
        do_set_parameter(key, params.at(key));
      }
    }
  }
};

} // namespace System
} // namespace ScriptInterface

#endif
