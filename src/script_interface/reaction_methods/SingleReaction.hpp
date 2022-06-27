/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_REACTION_METHODS_SINGLE_REACTION_HPP
#define SCRIPT_INTERFACE_REACTION_METHODS_SINGLE_REACTION_HPP

#include "core/reaction_methods/SingleReaction.hpp"
#include "script_interface/ScriptInterface.hpp"
#include <numeric>
#include <vector>

namespace ScriptInterface {
namespace ReactionMethods {

class SingleReaction : public AutoParameters<SingleReaction> {
public:
  SingleReaction() {
    add_parameters({
        {"gamma", AutoParameter::read_only, [this]() { return m_sr->gamma; }},
        {"reactant_types", AutoParameter::read_only,
         [this]() { return m_sr->reactant_types; }},
        {"reactant_coefficients", AutoParameter::read_only,
         [this]() { return m_sr->reactant_coefficients; }},
        {"product_types", AutoParameter::read_only,
         [this]() { return m_sr->product_types; }},
        {"product_coefficients", AutoParameter::read_only,
         [this]() { return m_sr->product_coefficients; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    m_sr = std::make_shared<::ReactionMethods::SingleReaction>(
        get_value<double>(params, "gamma"),
        get_value<std::vector<int>>(params, "reactant_types"),
        get_value<std::vector<int>>(params, "reactant_coefficients"),
        get_value<std::vector<int>>(params, "product_types"),
        get_value<std::vector<int>>(params, "product_coefficients"));
  }

  std::shared_ptr<::ReactionMethods::SingleReaction> get_reaction() {
    return m_sr;
  }

private:
  std::shared_ptr<::ReactionMethods::SingleReaction> m_sr;
};
} /* namespace ReactionMethods */
} /* namespace ScriptInterface */

#endif