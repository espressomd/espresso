/*
 * Copyright (C) 2022 The ESPResSo project
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

#include "bond_breakage.hpp"

#include "script_interface/ScriptInterface.hpp"

#include <memory>

namespace ScriptInterface {
namespace BondBreakage {

class BreakageSpec : public AutoParameters<BreakageSpec> {
public:
  BreakageSpec() : m_breakage_spec(new ::BondBreakage::BreakageSpec) {
    add_parameters({
        {"breakage_length", m_breakage_spec->breakage_length},
        {"action_type",
         [this](const Variant &v) {
           m_breakage_spec->action_type =
               ::BondBreakage::ActionType(boost::get<int>(v));
         },
         [this]() { return Variant(int(m_breakage_spec->action_type)); }},
    });
  }

  std::shared_ptr<::BondBreakage::BreakageSpec> breakage_spec() const {
    return m_breakage_spec;
  }

private:
  std::shared_ptr<::BondBreakage::BreakageSpec> m_breakage_spec;
};

} // namespace BondBreakage
} // namespace ScriptInterface
