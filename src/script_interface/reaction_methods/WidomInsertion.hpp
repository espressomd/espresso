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
  std::shared_ptr<::ReactionMethods::ReactionAlgorithm> const
  RE() const override {
    return m_re;
  }

  WidomInsertion() {
    add_parameters(
        {{"search_algorithm",
          [this](Variant const &) { throw_on_exclusion_change(); },
          []() { return none; }},
         {"exclusion_range",
          [this](Variant const &) { throw_on_exclusion_change(); },
          [this]() { return m_exclusion->get_parameter("exclusion_range"); }},
         {"exclusion_radius_per_type",
          [this](Variant const &) { throw_on_exclusion_change(); },
          [this]() {
            return m_exclusion->get_parameter("exclusion_radius_per_type");
          }}});
  }

  void do_construct(VariantMap const &params) override {
    setup_neighbor_search(params);
    context()->parallel_try_catch([&]() {
      auto exclusion = m_exclusion->get_instance();
      if (exclusion->exclusion_range != 0. or
          not exclusion->exclusion_radius_per_type.empty()) {
        throw std::runtime_error("No search algorithm for WidomInsertion");
      }
      m_re = std::make_shared<::ReactionMethods::WidomInsertion>(
          context()->get_comm(), get_value<int>(params, "seed"),
          get_value<double>(params, "kT"), m_exclusion->get_instance());
    });
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "calculate_particle_insertion_potential_energy") {
      Variant result;
      context()->parallel_try_catch([&]() {
        auto const reaction_id = get_value<int>(params, "reaction_id");
        auto const index = get_reaction_index(reaction_id);
        result = m_re->calculate_particle_insertion_potential_energy(index);
      });
      return result;
    }
    return ReactionAlgorithm::do_call_method(name, params);
  }

private:
  std::shared_ptr<::ReactionMethods::WidomInsertion> m_re;
  void throw_on_exclusion_change() const {
    if (context()->is_head_node()) {
      throw std::runtime_error("No search algorithm for WidomInsertion");
    }
  }
};

} /* namespace ReactionMethods */
} /* namespace ScriptInterface */

#endif
