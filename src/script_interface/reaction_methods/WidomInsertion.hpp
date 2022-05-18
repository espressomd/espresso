/*
 * Copyright (C) 2021 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_REACTION_METHODS_WIDOM_INSERTION_HPP
#define SCRIPT_INTERFACE_REACTION_METHODS_WIDOM_INSERTION_HPP

#include "ReactionAlgorithm.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/reaction_methods/ReactionAlgorithm.hpp"
#include "core/reaction_methods/WidomInsertion.hpp"

#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace ReactionMethods {

class WidomInsertion : public ReactionAlgorithm {
public:
  std::shared_ptr<::ReactionMethods::ReactionAlgorithm> RE() override {
    return m_re;
  }

  WidomInsertion() {
    add_parameters({{"search_algorithm",
                     [](Variant const &) {
                       throw std::runtime_error(
                           "No search algorithm for WidomInsertion");
                     },
                     []() { return none; }}});
  }

  void do_construct(VariantMap const &params) override {
    m_re = std::make_shared<::ReactionMethods::WidomInsertion>(
        get_value<int>(params, "seed"), get_value<double>(params, "kT"), 0.,
        std::unordered_map<int, double>{});
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &parameters) override {
    if (name == "calculate_particle_insertion_potential_energy") {
      auto const reaction_id = get_value<int>(parameters, "reaction_id");
      auto const index = get_reaction_index(reaction_id);
      auto &reaction = *m_reactions[index]->get_reaction();
      return m_re->calculate_particle_insertion_potential_energy(reaction);
    }
    return ReactionAlgorithm::do_call_method(name, parameters);
  }

private:
  std::shared_ptr<::ReactionMethods::WidomInsertion> m_re;
};
} /* namespace ReactionMethods */
} /* namespace ScriptInterface */

#endif