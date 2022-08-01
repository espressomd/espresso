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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_REACTION_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_REACTION_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include "EKReactant.hpp"
#include "LatticeWalberla.hpp"
#include "optional_reduction.hpp"

#include "walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp"
#include "walberla_bridge/electrokinetics/reactions/EKReactionImplBulk.hpp"
#include "walberla_bridge/electrokinetics/reactions/EKReactionIndexed.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameter.hpp"

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {
class EKReaction : public AutoParameters<::walberla::EKReactionBase> {
public:
  [[nodiscard]] virtual std::shared_ptr<::walberla::EKReactionBase>
  get_instance() const = 0;

  [[nodiscard]] Utils::Vector3i
  get_mapped_index(const Utils::Vector3i &node) const {
    auto output = node;
    const auto shape = get_instance()->get_lattice()->get_grid_dimensions();
    for (auto i : {0, 1, 2}) {
      if (node[i] < 0) {
        output[i] = node[i] + shape[i];
      }

      if (output[i] < 0 or output[i] >= shape[i]) {
        auto constexpr formatter = Utils::Vector3i::formatter(", ");
        std::stringstream ss;
        ss << "provided index [" << formatter << node
           << "] is out of range for shape [" << formatter << shape << "]\n";
        throw std::runtime_error(ss.str());
      }
    }

    return output;
  }
};

class EKBulkReaction : public EKReaction {
public:
  void do_construct(VariantMap const &args) override {

    auto lattice =
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice();

    auto reactant = get_value<std::vector<Variant>>(args, "reactants");
    std::vector<std::shared_ptr<::walberla::EKReactant>> output(
        reactant.size());
    auto get_instance = [](Variant &ekreactant) {
      return get_value<std::shared_ptr<EKReactant>>(ekreactant)->get_instance();
    };
    std::transform(reactant.begin(), reactant.end(), output.begin(),
                   get_instance);

    m_ekreaction = std::make_shared<::walberla::EKReactionImplBulk>(
        lattice, output, get_value<double>(args, "coefficient"));

    add_parameters({{"coefficient",
                     [this](Variant const &v) {
                       m_ekreaction->set_coefficient(get_value<double>(v));
                     },
                     [this]() { return m_ekreaction->get_coefficient(); }}});
  }

  [[nodiscard]] std::shared_ptr<::walberla::EKReactionBase>
  get_instance() const override {
    return m_ekreaction;
  }

private:
  std::shared_ptr<::walberla::EKReactionImplBulk> m_ekreaction;
};

class EKIndexedReaction : public EKReaction {
public:
  void do_construct(VariantMap const &args) override {

    auto lattice =
        get_value<std::shared_ptr<LatticeWalberla>>(args, "lattice")->lattice();

    auto reactant = get_value<std::vector<Variant>>(args, "reactants");
    std::vector<std::shared_ptr<::walberla::EKReactant>> output(
        reactant.size());
    auto get_instance = [](Variant &ekreactant) {
      return get_value<std::shared_ptr<EKReactant>>(ekreactant)->get_instance();
    };
    std::transform(reactant.begin(), reactant.end(), output.begin(),
                   get_instance);

    m_ekreaction = std::make_shared<::walberla::EKReactionIndexed>(
        lattice, output, get_value<double>(args, "coefficient"));

    add_parameters(
        {{"coefficient",
          [this](Variant const &v) {
            m_ekreaction->set_coefficient(get_value<double>(v));
          },
          [this]() { return m_ekreaction->get_coefficient(); }},
         {"shape", AutoParameter::read_only, [this]() {
            return m_ekreaction->get_lattice()->get_grid_dimensions();
          }}});
  }

  [[nodiscard]] Variant do_call_method(std::string const &method,
                                       VariantMap const &parameters) override {
    if (method == "set_node_is_boundary") {
      m_ekreaction->set_node_is_boundary(
          get_mapped_index(get_value<Utils::Vector3i>(parameters, "node")),
          get_value<bool>(parameters, "is_boundary"));
      return none;
    }
    if (method == "get_node_is_boundary") {
      auto const result = m_ekreaction->get_node_is_boundary(
          get_mapped_index(get_value<Utils::Vector3i>(parameters, "node")));
      return optional_reduction_with_conversion(result);
    }
    return none;
  }

  [[nodiscard]] std::shared_ptr<::walberla::EKReactionBase>
  get_instance() const override {
    return m_ekreaction;
  }

private:
  std::shared_ptr<::walberla::EKReactionIndexed> m_ekreaction;
};
} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
#endif
