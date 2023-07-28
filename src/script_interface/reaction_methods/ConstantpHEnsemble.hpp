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

#ifndef SCRIPT_INTERFACE_REACTION_METHODS_CONSTANT_PH_HPP
#define SCRIPT_INTERFACE_REACTION_METHODS_CONSTANT_PH_HPP

#include "ReactionAlgorithm.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/reaction_methods/ConstantpHEnsemble.hpp"
#include "core/reaction_methods/ReactionAlgorithm.hpp"
#include "core/reaction_methods/utils.hpp"

#include <memory>
#include <unordered_map>

namespace ScriptInterface {
namespace ReactionMethods {

class ConstantpHEnsemble : public ReactionAlgorithm {
public:
  std::shared_ptr<::ReactionMethods::ReactionAlgorithm> RE() override {
    return m_re;
  }
  std::shared_ptr<::ReactionMethods::ReactionAlgorithm> const
  RE() const override {
    return m_re;
  }

  ConstantpHEnsemble() {
    add_parameters({
        {"constant_pH",
         [this](Variant const &v) {
           m_re->m_constant_pH = get_value<double>(v);
         },
         [this]() { return m_re->m_constant_pH; }},
    });
  }

  void do_construct(VariantMap const &params) override {
    setup_neighbor_search(params);
    context()->parallel_try_catch([&]() {
      m_re = std::make_shared<::ReactionMethods::ConstantpHEnsemble>(
          context()->get_comm(), get_value<int>(params, "seed"),
          get_value<double>(params, "kT"), m_exclusion->get_instance(),
          get_value<double>(params, "constant_pH"));
    });
  }

protected:
  double calculate_factorial_expression(
      ::ReactionMethods::SingleReaction const &reaction,
      std::unordered_map<int, int> const &particle_numbers) const override {
    return ::ReactionMethods::calculate_factorial_expression_cpH(
        reaction, particle_numbers);
  }

private:
  std::shared_ptr<::ReactionMethods::ConstantpHEnsemble> m_re;
};
} /* namespace ReactionMethods */
} /* namespace ScriptInterface */

#endif
