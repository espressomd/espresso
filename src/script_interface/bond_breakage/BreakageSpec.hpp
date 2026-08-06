/*
 * Copyright (C) 2022-2026 The ESPResSo project
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
#include <stdexcept>
#include <string>
#include <variant>

namespace ScriptInterface {
namespace BondBreakage {

class BreakageSpec : public AutoParameters<BreakageSpec> {
public:
  BreakageSpec()
      : m_breakage_spec(std::make_shared<::BondBreakage::BreakageSpec>()) {
    add_parameters({
        {"breakage_length", m_breakage_spec->breakage_length},
        {"action_type",
         [this](Variant const &v) {
           m_breakage_spec->action_type = ::BondBreakage::ActionType{
               m_breakage_str_to_enum.at(std::get<std::string>(v))};
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
  void do_construct(VariantMap const &params) override {
    context()->parallel_try_catch([&]() {
      // all parameters are required: reject construction if any is omitted,
      // because "breakage_length" has no meaningful default (0.0 is a valid
      // value that would silently queue every bond of that type for breakage)
      for (auto const &key : valid_parameters()) {
        if (not params.contains(std::string{key})) {
          throw std::runtime_error("Parameter '" + std::string{key} +
                                   "' is missing");
        }
      }
    });
    AutoParameters<BreakageSpec>::do_construct(params);
  }

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
