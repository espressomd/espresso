/*
 * Copyright (C) 2026 The ESPResSo project
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

#include <config/config.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameter.hpp>
#include <script_interface/system/Leaf.hpp>

#include <memory>
#include <optional>
#include <string>

namespace ScriptInterface::EK {

class Container : public AutoParameters<Container, System::Leaf> {
public:
  Container() {
    add_parameters({
        {"solver",
         [](Variant const &v) {
           if (not is_none(v)) {
             throw WriteError("solver");
           }
         },
         []() { return Variant{None{}}; }},
    });
  }
};

} // namespace ScriptInterface::EK
