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

#include "core/bond_breakage/bond_breakage.hpp"

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
           m_breakage_spec->action_type = ::BondBreakage::ActionType{
               m_breakage_str_to_enum.at(boost::get<std::string>(v))};
         },
         [this]() {
           return Variant(
               m_breakage_enum_to_str.at(m_breakage_spec->action_type));
         }},
    });
  }

  std::shared_ptr<::BondBreakage::BreakageSpec> breakage_spec() const {
    return m_breakage_spec;
  }

private:
  std::shared_ptr<::BondBreakage::BreakageSpec> m_breakage_spec;
  std::unordered_map<::BondBreakage::ActionType, std::string>
      m_breakage_enum_to_str = {
          {::BondBreakage::ActionType::NONE, "none"},
          {::BondBreakage::ActionType::DELETE_BOND, "delete_bond"},
          {::BondBreakage::ActionType::REVERT_BIND_AT_POINT_OF_COLLISION,
           "revert_bind_at_point_of_collision"}};
  std::unordered_map<std::string, ::BondBreakage::ActionType>
      m_breakage_str_to_enum = {
          {"none", ::BondBreakage::ActionType::NONE},
          {"delete_bond", ::BondBreakage::ActionType::DELETE_BOND},
          {"revert_bind_at_point_of_collision",
           ::BondBreakage::ActionType::REVERT_BIND_AT_POINT_OF_COLLISION}};
};

} // namespace BondBreakage
} // namespace ScriptInterface
